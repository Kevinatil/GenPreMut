from dataclasses import dataclass
from typing import Any, Callable, Dict, List, NewType, Optional, Tuple, Union
from collections.abc import Mapping

import numpy as np
import torch

from transformers.data.data_collator import PreTrainedTokenizerBase


@dataclass
class DataCollatorForMaskedGeneration:

    tokenizer: PreTrainedTokenizerBase
    mlm_probability: Any = 0.15
    max_mask: int = 5
    device: Any = 'cuda'
    pad_to_multiple_of: Optional[int] = None
    tf_experimental_compile: bool = False
    return_tensors: str = "pt"


    def __call__(self, features):
        return self.torch_call(features)

    def torch_call(self, examples: List[Union[List[int], Any, Dict[str, Any]]]) -> Dict[str, Any]:
        batch = self.tokenizer.pad(examples, return_tensors="pt", pad_to_multiple_of=self.pad_to_multiple_of)

        special_tokens_mask = batch.pop("special_tokens_mask", None)
        input_ids, len_ = self.torch_mask_tokens(
            batch["input_ids"], special_tokens_mask=special_tokens_mask
        )
        batch["token_ids"] = batch["input_ids"] # no mask
        batch["input_ids"] = input_ids # has mask
        for key, value in batch.items():
            batch[key] = value[:len_].to(self.device)
        return batch

    def torch_mask_tokens(self, inputs: Any, special_tokens_mask: Optional[Any] = None) -> Tuple[Any, Any]:
        """
        Prepare masked tokens inputs/labels for masked language modeling: 80% MASK, 10% random, 10% original.
        """
        while True:
            if type(self.mlm_probability) == float:
                # We sample a few tokens in each sequence for MLM training (with probability `self.mlm_probability`)
                probability_matrix = torch.full(inputs.shape, self.mlm_probability)
            else:
                probability_matrix = self.mlm_probability.repeat(inputs.shape[0], 1)

            if special_tokens_mask is None:
                special_tokens_mask = [
                    self.tokenizer.get_special_tokens_mask(val, already_has_special_tokens=True) for val in inputs.tolist()
                ]
                special_tokens_mask = torch.tensor(special_tokens_mask, dtype=torch.bool)
            else:
                special_tokens_mask = special_tokens_mask.bool()

            probability_matrix.masked_fill_(special_tokens_mask, value=0.0)
            masked_indices = torch.bernoulli(probability_matrix).bool()
            # labels[~masked_indices] = -100  # We only compute loss on masked tokens

            inputs[masked_indices] = self.tokenizer.convert_tokens_to_ids(self.tokenizer.mask_token)

            masked_indices = (masked_indices.sum(dim=1) < self.max_mask) & (masked_indices.sum(dim=1) >= 1)

            if masked_indices.any():
                return inputs[masked_indices], masked_indices.sum().item()
