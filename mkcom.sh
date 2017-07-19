dir=`pwd`
newdir=`echo $dir | sed -e 's:/lab/work/arushiv/::g'`


if [ $# -eq 0 ]; then
	dirname="/home/arushiv/$newdir"
else
	dirname="${1}/$newdir"
fi
    
mkdir -p $dirname
cp /home/arushiv/toolScripts/commands.py $dirname
ln -s $dirname/commands.py .

# git init $dirname
# git add $dirname/commands.py

