#!/bin/sh
set -e +x

openssl aes-256-cbc \
  -K $encrypted_da75db262454_key \
  -iv $encrypted_da75db262454_iv \
  -in tests/data/cases/ytk_device.tar.xz.enc \
  -out tests/data/cases/ytk_device.tar.xz -d
