#!/usr/bin/env bash
# Linux/OSX script to run sphinx-apidoc with appropriate options
printf -- '-- Removing current apidoc content\n'
rm -rf apidoc
printf -- '-- Running sphinx-apidoc command\n'
export SPHINX_APIDOC_OPTIONS=members,show-inheritance
sphinx-apidoc ../proteuslib "../proteuslib/*tests" -o apidoc