#!/bin/bash
# Script runs commands in parallel
for task in "$@"; do {
  $task &
} done
