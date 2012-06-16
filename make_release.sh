#!/bin/bash

VER=`cat ../trunk/VERSION`
echo "Making rQuant release version $VER"
export RDIR=rQuant-$VER;
mkdir -p $RDIR
rsync -a ../trunk/ $RDIR
rm -rf $RDIR/tmp
rm -rf $RDIR/examples
svn up $RDIR/examples
rm -rf $RDIR/mosek
svn up $RDIR/mosek
rm -rf $RDIR/mex
svn up $RDIR/mex
rm -rf $RDIR/test_data*
svn up $RDIR/test_data
find $RDIR -name .svn -exec rm -rf {} \;
find $RDIR -name sge -exec rm -rf {} \;
find $RDIR -name \*~ -exec rm {} \;
rm -f $RDIR/octave-core $RDIR/*.mat
rm -f $RDIR/bin/rquant_config.sh
rm -f $RDIR/bin/rquant_config.sh.bak
sed -e s/"RQUANT_VERSION=X.X"/"RQUANT_VERSION=$VER"/g $RDIR/bin/rquant_config.sh.template  > tmp_file
chmod 755 tmp_file
mv tmp_file $RDIR/bin/rquant_config.sh.template
mv $RDIR/bin/rquant_config.sh.template $RDIR/bin/rquant_config.sh
tar -cjvf $RDIR.tar.bz2 $RDIR
