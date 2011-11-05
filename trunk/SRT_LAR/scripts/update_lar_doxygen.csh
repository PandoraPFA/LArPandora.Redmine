#!/bin/csh

setenv LARCODE /grid/fermiapp/lbne/lar/code/

setenv LAR_WWW /nusoft/app/web/htdoc/larsoft/doxsvn

echo Running Doxygen
cd $LARCODE
rm -rf doxygen/dox/html
mkdir doxygen/dox/html
./doxygen/doxygen-1.7.1/bin/doxygen doxygen/doxylar >! lar_doxygen.log

cd $LAR_WWW
rm -rf ./html
echo Copying output
cp -r $LARCODE/doxygen/dox/html/ ./html

echo done at `date`


