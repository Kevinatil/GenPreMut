import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM, Trainer, TrainingArguments, DataCollatorForLanguageModeling
from datasets import Dataset

save_folder = "ckpt/finetune"
batch_size = 16

init_data_path = ""


def get_training_sequences(path):
    seqs = []
    f = open(path, 'r')
    for line in f:
        seqs.append(line.strip())

    return seqs



if __name__ == "__main__":
    tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
    model = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t33_650M_UR50D")
    data_collator = DataCollatorForLanguageModeling(tokenizer=tokenizer, mlm_probability=0.15)


    train_sequences = get_training_sequences(init_data_path)
    train_tokenized = tokenizer(train_sequences)
    train_dataset = Dataset.from_dict(train_tokenized)
    train_args = TrainingArguments(
        output_dir=save_folder,
        save_strategy = "epoch",
        learning_rate=1e-4,
        per_device_train_batch_size=batch_size,
        per_device_eval_batch_size=batch_size,
        num_train_epochs=5,
        weight_decay=0.01,
        warmup_steps=1000,
        report_to="none",
    )


    trainer = Trainer(
        model,
        train_args,
        train_dataset=train_dataset,
        tokenizer=tokenizer,
        data_collator=data_collator,
    )
    trainer.train()