#!/usr/bin/env bash

if ! [ -f hg38/hg38.fa.sizes ]; then
        genomepy install -p UCSC -g ./ -r "^chr[\dXYM]{1,2}$" -f hg38 -m hard
fi

if ! [ -f mm10/mm10.fa.sizes ]; then
        genomepy install -p UCSC -g ./ -r "^chr[\dXYM]{1,2}$" -f mm10 -m hard
fi