#!/usr/bin/env bash
#Simple installation script for MuTools. Simply appends aliases to each
#tool at the end of the user's home directory.

echo '#The following are aliases to each MuTools utility.' >> ~/.bashrc
echo "alias multimutect=$(pwd)/multimutect/multimutect.py" >> ~/.bashrc
for i in $( echo premutect/*.py ); do
   base=$(basename $i .py)
   if [ "$base" = "__init__" ]
      then continue
   else
      echo "alias $base=$(pwd)/$i" >> ~/.bashrc
   fi
done
echo "alias combine=$(pwd)/postmutect/combine.py" >> ~/.bashrc
echo "alias catenate=$(pwd)/postmutect/catenate.py" >> ~/.bashrc
