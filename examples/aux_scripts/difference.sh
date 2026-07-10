#!/bin/bash

DIFFTOOL=meld

$DIFFTOOL output/*.rmta.out reference/*.rmta.ref
