#!/bin/sh
set -e +x

for file in tests/data/*.enc; do
  echo "Decrypting $file..."
  openssl aes-256-cbc              \
    -K $encrypted_3f3dad14a293_key \
    -iv $encrypted_3f3dad14a293_iv \
    -in "$file"                    \
    -out "${file%.enc}"           \
    -d
done
