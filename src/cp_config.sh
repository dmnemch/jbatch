#!/bin/bash

SRC_CONFIG_FILE="$(dirname "$0")/../config/jb_config.yaml"
DEST_CONFIG_DIR="$HOME/.config/jb"
DEST_CONFIG_FILE="$DEST_CONFIG_DIR/jb_config.yaml"

if [ ! -d "$DEST_CONFIG_DIR" ]; then
  mkdir -p "$DEST_CONFIG_DIR"
fi

cp "$SRC_CONFIG_FILE" "$DEST_CONFIG_FILE"
echo "\033[32mGenerated default config in $DEST_CONFIG_DIR, \033[31mbut it lacks crucial argument (e.g. partition) \033[0m -- please add them manually"
