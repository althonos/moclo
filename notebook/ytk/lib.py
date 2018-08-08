# coding: utf-8
import collections
import datetime
import re
import textwrap

import fs
import Bio.SeqIO
import Bio.SeqUtils
import pandas
import qgrid
import six
from notebook import notebookapp
from IPython.display import display, HTML

from moclo.kits import ytk
from moclo.registry.base import CombinedRegistry
from moclo.registry.ytk import YTKRegistry, PTKRegistry


def load_registry(ptk=False):
    registry = CombinedRegistry() << YTKRegistry()
    if ptk:
        registry << PTKRegistry()
    return registry

def ask_plasmids():
    print(textwrap.dedent("""
    Fill the table with as many plasmids as desired, using the "Add Row" button
    to add more rows. When finished, run the next cell.
    """))
    df_plasmids = pandas.DataFrame({'id': ['psXXX'], 'name': ['plasmid_name']})
    qgrid_plasmids = qgrid.QgridWidget(df=df_plasmids, show_toolbar=True)
    display(qgrid_plasmids)
    return qgrid_plasmids

def ask_parts(registry, qgrid_plasmids):
    print(textwrap.dedent("""
    Use the selector to choose which parts to use in each plasmid. Parts are
    sorted by YTK type. When finished, run the next cell.
    """))
    # Extract user plasmids IDs and names
    df_plasmids = qgrid_plasmids.get_changed_df()
    names = {row['id']: row['name'] for _, row in df_plasmids.iterrows()}
    # Input parts dataframe
    types = {
        re.search('YTKPart(.*)', cls.__name__).group(1): cls
        for cls in ytk.YTKPart.__subclasses__()
    }
    # Create a small selected with part ID and name
    # for each possible part
    parts = {}
    for colname, part_type in types.items():
        categories = {
            '{} - {}'.format(p.entity.record.id, p.entity.record.description)
            for p in registry.values()
            if isinstance(p.entity, part_type)
        }
        parts[colname] = pandas.Categorical(
            [''] * len(names),
            categories=[''] + sorted(categories)
        )
    # Create the dataframe
    df_parts = pandas.DataFrame(parts)
    df_parts['id'] = df_plasmids['id']
    df_parts['name'] = df_plasmids['name']
    qgrid_parts = qgrid.QgridWidget(df=df_parts, show_toolbar=False)
    display(qgrid_parts)
    return qgrid_parts


def validate_assemblies(registry, qgrid_parts):
    print(textwrap.dedent("""
    Making sure all assemblies are valid:
    """))
    full = []
    # Extract
    for _, row in qgrid_parts.get_changed_df().sort_index().iterrows():
        *parts, plasmid_id, plasmid_name = filter(None, row)
        for part_id, part_name in (p.split(' - ') for p in parts):
            part = registry[part_id]
            idx = int(re.search('pYTK(\d\d\d)', part_id).group(1)) - 1
            full.append({
                'Plasmid Name': plasmid_name,
                'Plasmid ID': plasmid_id,
                'Part ID': part_id,
                'Part Position': "{}{}".format(chr(ord('A') + idx // 12), idx % 12 + 1),
                'Part Name': part.entity.record.description,
                'Resistance': part.resistance,
                'Record': part.entity.record,
                'Typed part': part.entity,
            })

    assemblies = collections.defaultdict(list)
    for p in full:
        assemblies[p['Plasmid ID']].append(p['Typed part'])
    for pid, parts in assemblies.items():
        print("- Validating plasmid", pid, end="... ")
        vector = parts.pop(next(i for i in range(len(parts)) if isinstance(parts[i], ytk.YTKCassetteVector)))
        assembly = vector.assemble(*parts)
        print("✓")

    df_full = pandas.DataFrame(full)
    df_full = df_full[['Plasmid ID', 'Plasmid Name', 'Part ID', 'Part Name', 'Resistance', 'Record', 'Typed part']]

    display(df_full
                .drop(columns=['Record', 'Typed part'])
                .sort_values("Plasmid ID")
                .reset_index()
                .drop(columns='index')
    )

    return df_full


def ask_concentrations(assemblies):
    print(textwrap.dedent("""
    Use the table below to input mass concentration for all your parts.
    """))


    df_strains = assemblies.drop(columns=['Plasmid Name', 'Plasmid ID'])\
                    .drop_duplicates(subset='Part ID')\
                    .sort_values('Part ID')\
                    .reset_index()\
                    .drop(columns='index')

    def part_size(part):
        return len(part.target_sequence())

    df_strains['Plasmid size'] = df_strains['Record'].apply(len)
    df_strains['Part size'] = df_strains['Typed part'].apply(part_size)
    df_strains['Sample concentration (ng/µL)'] = [0] * len(df_strains['Part ID'])

    concentration_widget = qgrid.QGridWidget(df=df_strains[['Part ID', 'Sample concentration (ng/µL)']])
    display(concentration_widget)
    return df_strains, concentration_widget

def dilution_table(df_strains, concentration_widget):
    df_strains['Sample concentration (ng/µL)'] = concentration_widget.get_changed_df().sort_index()['Sample concentration (ng/µL)']
    df_strains['MW (Dalton)'] = df_strains['Record'].apply(lambda r: Bio.SeqUtils.molecular_weight(r.seq, 'DNA', True, True))
    df_strains['Sample concentration (fmol/µL)'] = 10 ** 6 * df_strains['Sample concentration (ng/µL)'] / df_strains['MW (Dalton)']
    df_strains['µL of sample to 40fmol'] = 20 / df_strains['Sample concentration (fmol/µL)']
    df_strains['Dilution factor to 20fmol/µL'] = round(1 / df_strains['µL of sample to 40fmol']).apply(int)
    return df_strains.drop(columns=['Record', 'Typed part'])


def generate_gb_files(assemblies):


    now = datetime.datetime.now()

    filename = "MoClo - {}.zip".format(datetime.datetime.now())
    with fs.open_fs('zip://{}'.format(filename), create=True) as archive:

        for plasmid_id in set(assemblies['Plasmid ID']):

            vector = None
            modules = []

            for _, row in assemblies.iterrows():
                if row['Plasmid ID'] == plasmid_id:
                    if isinstance(row['Typed part'], ytk.YTKCassetteVector):
                        vector = row['Typed part']
                    else:
                        modules.append(row['Typed part'])

            record = vector.assemble(*modules)
            record.name = row['Plasmid Name']
            record.id = row['Plasmid ID']

            with archive.open('{}.gb'.format(plasmid_id), 'w') as f:
                Bio.SeqIO.write(record, f, 'gb')


    server = next(notebookapp.list_running_servers())
    display(HTML("<a href='{0}files/{1}'>{1}<a>".format(server['url'], filename)))
