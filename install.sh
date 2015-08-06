#!/usr/bin/env bash
#Simple symlink installation script for MuTools.


ln -s $(pwd)/multimutect/multimutect.py /usr/bin/multimutect
for i in $( echo premutect/*.py ); do
   if [ "$(basename $i)" = '__init__.py' ]
      then continue
   else
      ln -s $(pwd)/$i /usr/bin/$(basename -s .py  $i)
   fi
done
ln -s $(pwd)/postmutect/combine.py /usr/bin/combine
ln -s $(pwd)/postmutect/catenate.py /usr/bin/catenate
