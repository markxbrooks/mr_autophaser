# mr_autophaser
python script to solve a crystal structure by molecular replacement using Phaser

example usage:

```
ccp4-python ./mr_autophaser.py -m test/beta_blip.mtz -1 test/beta.pdb -2 test/blip.pdb
```

Unit test

```
ccp4-python -m pytest
```
