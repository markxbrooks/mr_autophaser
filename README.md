# mr_autophaser
python script to solve a crystal structure by molecular replacement using Phaser

example usage:

```
ccp4-python ./mr_autophaser.py -m test/beta_blip.mtz -1 test/beta.pdb -2 test/blip.pdb
```


```
optional arguments:
  -h, --help            show this help message and exit
  -m mtzin, --mtz_input mtzin
                        MTZ input file
  -1 PDBIN1, --pdb_input1 PDBIN1
                        pdb input file
  -2 PDBIN2, --pdb_input2 PDBIN2
                        pdb input file
  -3 PDBIN3, --pdb_input3 PDBIN3
                        pdb input file3
  -c chain, --chain chain
                        chain name, e.g. 'A','C',' ' Default is 'None'
  -n NUMBER1, --num_pdb1 NUMBER1
                        -num_pdb1, default is 1
  -o NUMBER2, --num_pdb2 NUMBER2
                        -num_pdb2, default is 1
  -p NUMBER3, --num_pdb3 NUMBER3
                        -num_pdb3, default is 1
```

```

Unit test

```
ccp4-python -m pytest
```
