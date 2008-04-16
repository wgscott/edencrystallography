#!/bin/csh

if ( $#argv < 2 ) then
    echo "Usage: mkeden.csh bindir executable
    exit
endif

cd ..
cat <<EOF > $1/eden
#!/bin/csh -f
setenv EDENHOME `pwd`
$2 \$*
EOF
chmod +x $1/eden

cat <<EOF > $1/ieden
#!/bin/csh -f
setenv EDENHOME `pwd`
\$EDENHOME/python/eden.py
EOF
chmod +x $1/ieden

