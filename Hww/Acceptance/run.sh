#!/bin/bash

INPUT=powhegzee.input
root -l -q plotZeeGen.C+\(\"$INPUT\"\)
root -l -q plotZeeReco.C+\(\"$INPUT\"\)
root -l -q plotZeeSelect.C+\(\"$INPUT\"\)

INPUT=fsrzee.input
#root -l -q plotZeeGen.C+\(\"$INPUT\"\)
#root -l -q plotZeeReco.C+\(\"$INPUT\"\)
#root -l -q plotZeeSelect.C+\(\"$INPUT\"\)

INPUT=bornzeeval.input
#root -l -q plotZeeGen.C+\(\"$INPUT\"\)
#root -l -q plotZeeReco.C+\(\"$INPUT\"\)
#root -l -q plotZeeSelect.C+\(\"$INPUT\"\)

INPUT=powhegzeeval.input
#root -l -q plotZeeGen.C+\(\"$INPUT\"\)
#root -l -q plotZeeReco.C+\(\"$INPUT\"\)
#root -l -q plotZeeSelect.C+\(\"$INPUT\"\)

INPUT=fsrzeeval.input
#root -l -q plotZeeGen.C+\(\"$INPUT\"\)
#root -l -q plotZeeReco.C+\(\"$INPUT\"\)
#root -l -q plotZeeSelect.C+\(\"$INPUT\"\)

INPUT=pythiapuzee.input
#root -l -q plotZeeGen.C+\(\"$INPUT\"\)
#root -l -q plotZeeReco.C+\(\"$INPUT\"\)
#root -l -q plotZeeSelect.C+\(\"$INPUT\"\)

INPUT=ewkzee.input
#root -l -q plotZeeGen.C+\(\"$INPUT\"\)
#root -l -q plotZeeReco.C+\(\"$INPUT\"\)
#root -l -q plotZeeSelect.C+\(\"$INPUT\"\)

INPUT=ewkzeeval.input
#root -l -q plotZeeGen.C+\(\"$INPUT\"\)
#root -l -q plotZeeReco.C+\(\"$INPUT\"\)
#root -l -q plotZeeSelect.C+\(\"$INPUT\"\)
