for ((i=1;i<=60;i++))
do
{
    cd src$i/
    rm -rf input/
    rm -rf bin/
    rm -rf checkpoint/
    rm -rf parfile/
    rm -rf pyplot/
    rm *.dat
    rm *.log
    rm *.mplstyle
    rm *.sh
    rm *.py
}