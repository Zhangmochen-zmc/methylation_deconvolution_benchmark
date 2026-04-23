## Workflow

The execution is divided into 2 main steps. Please follow them in order. This pipeline requires the use of `methylbert`.

### Step 1: Fine-tuning

Fine-tune the pretrained MethylBERT model using the training and test datasets:

```bash
methylbert finetune \
    -c input/training_data/train_seq.kmers_fixed.csv \
    -t input/training_data/test_seq.kmers_fixed.csv \
    -o output/finetune_model/ \
    -p output/hg38_output/ \
    -s 160 \
    --loss focal_bce \
    -b 600 \
    -e 1000 \
    --lr 4e-4 \
    --with_cuda
```

### Step 2: Deconvolution

Run batch deconvolution using the provided script:

```bash
python MethylBERT_deconv.py <input_folder> <output_base_folder>
```

---

### Notes

*   **More Information**: [https://github.com/CompEpigen/methylbert.git](https://github.com/CompEpigen/methylbert.git)
