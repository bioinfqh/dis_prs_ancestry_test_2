
xgmix_dir=$1

for i in {1..22}
do
if [[ -f $xgmix_dir/chm_$i.pkl.gz ]]
then
else
wget -O $xgmix_dir/chm_$i.pkl.gz https://github.com/AI-sandbox/XGMix-models/blob/master/build37/missing_0/chm_$i.pkl.gz
fi
done
