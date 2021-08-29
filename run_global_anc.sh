customer_id=$1
input_file=$2

bash /scripts/get_ancestry.sh $input_file /scripts/ancestry_map_$customer_id.html /scripts/output_$customer_id.Q
python3 /scripts/read_ancestry_file.py /scripts/output_$customer_id.Q /scripts/ancestry_$customer_id.json


