rm -rf data
rm -rf /staging/hhao9/temposub_highdim_data.tar.gz
rm true_order_and_stages.json
python3 gen.py

tar -czf data.tar.gz data
mv /home/hhao9/temposub_highdim/data.tar.gz /staging/hhao9/temposub_highdim_data.tar.gz
# rm -rf data

# bash run.sh

# ls -lh /staging/hhao9/temposub_highdim_data.tar.gz
