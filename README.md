# S-HOT: Scalable High-Order Tucker Decomposition (WSDM 2017)
Jinoh Oh, Kijung Shin, Evangelos E. Papalexakis, Christos Faloutsos, and Hwanjo Yu

## Webpage
Visit our official [webpage](http://dm.postech.ac.kr/shot/)!

## Citing S-HOT
We encourage you to cite our paper if you have used our implementation in your work. You can use the following BibTex citation:
>
<dl><pre>
@inproceedings{Oh:2017:SHOT,
  author = {Jinoh Oh and Kijung Shin and Evangelos E. Papalexakis and Christos Faloutsos and Hwanjo Yu},
  title = {S-HOT: Scalable High-Order Tucker Decomposition},
  booktitle = {Proceedings of the 10th ACM International Conference on Web Search and Data Mining}, 
  series = {WSDM '17},
  year = {2017},
  location = {Cambridge, UK},
  doi = {10.1145/3018661.3018721},
  publisher = {ACM},
  address = {New York, NY, USA}
}
</pre></dl>

## Requirements
- C++
- gFortran
- Supports AVX/SSE instructions ([web page](https://software.intel.com/sites/landingpage/IntrinsicsGuide/))
- (Optional) The followings are used for evaluation, but are not mandatory to run S-HOT.
  * Python 2.7
  * Matplotlib
  * MATLAB is used for evaluation 
  * Tensor Toolbox 2.6 ([web page](http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html))

## How to run
### Input file format
The input file should be a list of non-zero entries where a line corresponds to a non-zero entry and attribute values for the non-zero entry are white-space separated. For example,
> `<attribute_1> <attribute_2> <attribute_3> ... <attribute_n> <target_value>`

Note that each attribute value should be integer, but `<target_value>` can be float. The following is a typical example for a 4-order tensor.

>
<dl><pre>
5395 9458 4819 5770 0.773400
9769 4765 3811 5978 0.626600
3460 8639 2865 6085 0.776800
415 4501 5780 5264 0.697500
6641 357 1049 2985 0.708100
8061 6456 3453 8593 0.675900
1891 2662 374 6067 0.514700
5013 6668 6522 1579 0.453100
4602 6812 1189 5062 0.998600
35 2957 9111 5609 0.423300
9233 9168 8991 9018 0.471200
4444 8895 634 4166 0.399400
7059 5334 7542 3928 0.405200
2663 9732 7561 6378 0.545900
4011 8190 3314 3603 0.997900
5638 4177 5991 7627 0.762700
9102 5875 5914 637 0.912900
710 2522 3458 3463 0.149700
8234 9134 4184 7925 0.086500
9866 9009 9983 2927 0.787400</pre></dl>

### Test your own tensor
Run `<install_path>/bin/Shot --size-mode <int> --size-rank <string> --iteration <int> <tensor_path> <model_path>` in bash shell. The detailed configurations are presented in the `Configuration` section.

### Reproduce the reported results
1. Type `cd <installpath>/reproduce` in bash shell.
2. Run `sudo python run.py` in bash shell.
3. Check `<installpath>/reproduce/fig/`

* NOTE1: Make sure MATLAB can be run on your equipment by simply typing `matlab`.
* NOTE2: Since the script use randomly-generated tensor, the results could be "slightly" different.
* NOTE3: In our machine setting, reproducing whole results takes 1-2 days.
* WARNING: I strongly encourage to turn-off your swap space by `sudo swapoff -a`. Allowing swap space could slightly enhance scalability of MATLAB-based implementations but definitely takes enourmous time while does not change the overall trend.

## Configuration 
1. `--plain`: run plain S-HOT implementation.
2. `--scan`: run scan-optimized S-HOT implementation (Default)
3. `--size-mode <int>`: the size of an input tensor.
4. `--size-rank <string>`: the size of rank of each mode after decomposition. `<string>` is formatted as `<int>:<int>:<int>...`, which is the list of integers with delimiter `:`. For example, `--size-rank 4:8:12:16` means the rank of n-mode after decomposition is 4, 8, 12, 16, respectively.
5. `--size-thread <int>`: the size of parallel threads.
6. `--use-modelloss`: use factor model changes as loss (Default)
7. `--use-coreloss`: use core loss
8. `--use-iteration`: use number of iterations as terminal condition (Default)
9. `--use-epsilon`: use loss change as terminal condition
10. `--iteration`: set the number of total iteration (Default = 1)
11. `--epsilon`: set the loss tolerance (Default = 0.01)
