#!/bin/sh

set -e +x


cd /py

echo ""
echo "> before_install"
pip install -U -r ci/requirements.txt

echo ""
echo "> install"
pip install ./moclo*

echo ""
echo "> before_script"
pip install -U -r tests/requirements.txt

echo ""
echo "> script"
green
