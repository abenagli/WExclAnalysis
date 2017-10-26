#!/bin/bash

function indent ()
{
    emacs -batch ${1} --eval '(setq-default indent-tabs-mode nil)' \
          --eval '(c-set-style "linux")' \
          --eval '(setq tab-width 4)' \
          --eval '(setq c-basic-offset 4)' \
          --eval '(indent-region (point-min) (point-max))' \
          -f save-buffer
}

if [ -z ${1} ]; then
    for file in `find ./*/*.h; find ./*/*.cc`;
    do
        indent $file
    done
else
    indent ${1}
fi


