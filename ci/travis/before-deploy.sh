#!/bin/sh

set -e +x

case "$(echo $TRAVIS_TAG | cut -d/ -f1)" in
  "ytk")   DIR=moclo-ytk;;
  "cidar") DIR=moclo-cidar;;
  *)       DIR=moclo;;
esac

cd $DIR
echo "Deploying $DIR..."

python setup.py check -rms sdist -d ../dist bdist_wheel -d ../dist

cd $TRAVIS_BUILD_DIR
