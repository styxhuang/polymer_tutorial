#!bin/bash
t=(100 150 200 300 400 500 600 700 800 900)
name=pressure
#group='1 2 3 4'
#group='5 6 7 9'
groupX='33'
groupY='37'
groupZ='41'
xgroup='17'
ygroup='18'
zgroup='19'

rm -r youngs_data-$name
mkdir youngs_data-$name

for i in ${t[@]}
do
    mkdir youngs_data-$name/T$i
    cd T$i
    for ii in $(ls -d sim*)
    do
	cd $ii
	pwd
	gmx energy -f y-x.edr -o $i-$ii-x-$name-frc.xvg -sum <<< $groupX
	gmx energy -f y-x.edr -o $i-$ii-x-$name-box.xvg <<< $xgroup
	sed -i '/^#/d' $i-$ii-x-$name-frc.xvg
	sed -i '/^@/d' $i-$ii-x-$name-frc.xvg
	sed -i '/^#/d' $i-$ii-x-$name-box.xvg
	sed -i '/^@/d' $i-$ii-x-$name-box.xvg
	paste $i-$ii-x-$name-box.xvg $i-$ii-x-$name-frc.xvg | awk '{print $1, $2, $4}' > $i-$ii-x-$name.xvg
	
	gmx energy -f y-y.edr -o $i-$ii-y-$name-frc.xvg -sum <<< $groupY
	gmx energy -f y-y.edr -o $i-$ii-y-$name-box.xvg <<< $ygroup
        sed -i '/^#/d' $i-$ii-y-$name-frc.xvg
        sed -i '/^@/d' $i-$ii-y-$name-frc.xvg
        sed -i '/^#/d' $i-$ii-y-$name-box.xvg
        sed -i '/^@/d' $i-$ii-y-$name-box.xvg
	paste $i-$ii-y-$name-box.xvg $i-$ii-y-$name-frc.xvg | awk '{print $1, $2, $4}' > $i-$ii-y-$name.xvg

	gmx energy -f y-z.edr -o $i-$ii-z-$name-frc.xvg -sum <<< $groupZ
	gmx energy -f y-z.edr -o $i-$ii-z-$name-box.xvg <<< $zgroup
        sed -i '/^#/d' $i-$ii-z-$name-frc.xvg
        sed -i '/^@/d' $i-$ii-z-$name-frc.xvg
        sed -i '/^#/d' $i-$ii-z-$name-box.xvg
        sed -i '/^@/d' $i-$ii-z-$name-box.xvg
        paste $i-$ii-z-$name-box.xvg $i-$ii-z-$name-frc.xvg | awk '{print $1, $2, $4}' > $i-$ii-z-$name.xvg

	for iii in $(ls *$name.xvg)
	do
	    mv $iii ../../youngs_data-$name/T$i/
	done
	cd ..
    done
    cd ..
    echo Finish
    pwd
done
cp -r youngs_data-$name $DROPBOX/0-dump/VE-Sty-T-effect/
