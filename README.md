﻿## 1.程序功能简介
本软件为GNSS接收机坐标估算软件，通过已给的导航定位源文件，可借助该程序进行导入数据并计算得到经过最小二乘法、卡尔曼滤波法、卡尔曼滤波CV模型与卡尔曼滤波CA模型四种算法改正后的GNSS坐标，并将结果进行绘图。软件的主要功能包括数据的导入与存储、解算算法的选择、解算历元的范围选择，解算状态的显示和解算结果的绘图。

### 1.1数据载入
文件类型：TXT格式；

数据输入格式为：头文件+数据文件 数据文件包括单次历元内的基本信息与历元内每颗卫星的卫星编号+伪距观测值+卫星钟差+卫星对流层改正数+卫星的坐标信息+卫星的高度角；

首先读取文件的头文件，提取测站的近似坐标，再读取个每个历元的基本信息，再通过基本信息循环读入每颗卫星的信息；

### 1.2坐标值的计算
程序先通过头文件确定接收机的近似坐标值，再通过卫星伪距单点定位方程列出伪距观测方程，并对伪距进行泰勒级数展开，之后通过相应的算法计算坐标的改正数，并将改正数与近似坐标相加，从而得到接收机坐标较为准确的估计值。

### 1.3数据信息的显示与绘图
在导入数据后，程序会自动提取头文件内容中的近似坐标，并显示在UI上，点击解算后，UI会实时显示正在解算第多少个历元，并将当前历元计算得到的改正数和改正后的坐标实时显示在UI上，同时将这些信息存储在矩阵变量内，在点击绘图后，会自动以历元为横坐标，解算结果为纵坐标进行绘图。

### 1.4数据的保存
点击解算后，UI会实时显示正在解算第多少个历元，并将当前历元计算得到的改正数和改正后的坐标实时显示在UI上，同时将这些信息存储在矩阵变量内，在点击保存后，选择要保存的路径和要保存的文件名，程序会自动以历元数+X坐标+Y坐标+Z坐标+X改正数+Y改正数+Z改正数的格式将矩阵中的内容输出到文件中。

## 2.主要算法设计与流程图
### 2.1算法设计
此程序是通过对卫星伪距单点定位方程进行线性化后，分别通过最小二乘算法、卡尔曼滤波算法、卡尔曼滤波CV算法和卡尔曼滤波CA算法计算坐标改正数后对接收机近似坐标进行改正，得到估计后的接收机坐标。

## 3.主要函数和变量说明
### 3.1主要函数
包括了读取数据、数据计算、数据绘图与数据输出的主要使用函数。

1. void **txt\_reader**(std::string filename, std::vector<Satellite>& vec)

功能：将TXT的数据读取到结构体数组中；

输入：数据文件的绝对路径，结构体变量名；

输出：接收机头文件中的近似坐标与所有历元的所有卫星信息

2. void MainWindow::**on\_LSM\_Button\_clicked**()

功能：运用最小二乘法计算接收机坐标的改正数；

输入：存储卫星信息的结构体数组与接收机近似坐标的数组

输出：通过最小二乘法解算的接收机坐标改正数与经过改正后的接收机坐标

3. void MainWindow::**on\_KALMAN\_Button\_clicked**()

功能：运用卡尔曼滤波法计算接收机坐标的改正数；

输入：存储卫星信息的结构体数组与接收机近似坐标的数组

输出：通过卡尔曼滤波法解算的接收机坐标改正数与经过改正后的接收机坐标

4. void MainWindow::**on\_KALMANCVButton\_clicked**()

功能：运用卡尔曼滤波CV模型计算接收机坐标的改正数；

输入：存储卫星信息的结构体数组与接收机近似坐标的数组

输出：通过卡尔曼滤波CV模型解算的接收机坐标改正数与改正后的接收机坐标

5. void MainWindow::**on\_KALMANCAButton\_clicked**()

功能：运用卡尔曼滤波CA模型计算接收机坐标的改正数；

输入：存储卫星信息的结构体数组与接收机近似坐标的数组

输出：通过卡尔曼滤波CA模型解算的接收机坐标改正数与改正后的接收机坐标

6. drawxx::**drawxx**(Eigen::MatrixXd POS1, int geshu1, double \*APPROX1, QWidget \*parent)

功能：对解算改正后的X坐标进行绘图

输入：经过改正后的坐标矩阵，解算的历元个数，接收机头文件中近似坐标

输出：X坐标的曲线图

7. drawyy::**drawyy**(Eigen::MatrixXd POS1, int geshu1, double \*APPROX1, QWidget \*parent)

功能：对解算改正后的Y坐标进行绘图

输入：经过改正后的坐标矩阵，解算的历元个数，接收机头文件中近似坐标

输出：Y坐标的曲线图

8. drawzz::**drawzz**(Eigen::MatrixXd POS1, int geshu1, double \*APPROX1, QWidget \*parent)

功能：对解算改正后的Z坐标进行绘图

输入：经过改正后的坐标矩阵，解算的历元个数，接收机头文件中近似坐标

输出：Z坐标的曲线图

9. void MainWindow::**baocun**()

功能：对解算结果进行保存

输入：改正数矩阵，经过改正后的坐标矩阵

输出：历元数+X坐标+Y坐标+Z坐标+X改正数+Y改正数+Z改正数的结果文件
### 3.2主要变量
1.typedef struct Data

{

`    `std::string PRN;

`    `double satnum, sats, satPos[3], sat\_Clock, Elevation, Azimuth, Prange[2], L[2], Trop\_Delay, Trop\_Map, Relativity, Sagnac, Tide\_Effect, Antenna\_Height, Sat\_Antenna, Tide\_Offset[2], Windup; 

}Data;//小结构体，用于存放单颗卫星的数据

2.typedef struct Satellite {

`    `int Satellite\_Num;

`    `double gpst, ztd, gps\_clock, glonass\_clock;

`    `std::vector<Data> data;

}Satellite;//大结构体，用于存放整个历元的数据

3.double APPROX\_POSITION[3];//全局变量，接收机近似坐标

1. double X0=0, Y0=0, Z0=0;//全局变量存储接收机位置

1. MatrixXd POS(2, 2);//坐标结果矩阵

1. MatrixXd ANS(4, 1); //坐标改正数矩阵

7.int geshu; //解算历元数

8.MatrixXd H0(1, 4);//创建系数矩阵

9.MatrixXd L(2, 2);//创建观测值矩阵

\10. MatrixXd lsm\_ans(4, 1);//创建改正数计算矩阵

\11. MatrixXd P0(2, 2);//创建权矩阵

\12. MatrixXd fuzhu(2, 2);//创建辅助计算矩阵

13.double dx; //创建X差值

14.double dy; //创建Y差值

15.double dz; //创建Z差值

17.double R0; //创建距离

18.double l; //创建系数l

19.double m; //创建系数m

20.double n; //创建系数n

## 4.程序界面

![github](https://github.com/Wang-Jie-Lucid-Sheep/Optimal_Estimation_C-/blob/main/Picture/Picture.png)
## 5.程序数据格式说明
### 5.1原始数据格式
数据文件格式：头文件+逐次历元的卫星数据

例如：

MARKER\_NAME：CUT0

INTERVAL： 30second

ANT\_TYPE：TRM59800.00     SCIS

REC\_TYPE：       TRIMBLE NETR9

APPROX\_POSITION：  -2364337.3850,   4870285.6069,  -3360809.7279  (m)

ANTENNA\_HEIGHT：         0.0000  (m)

Tropospheric Delay           : GPT Model

Tropospheric Mapping Function: Global Mapping Function

Navigation type:                 igs

Clock type: clock 30s

PRN,s: Satposition(X), Satposition(Y), Satposition(Z),   Sat Clock(m),   Elevation(°),     Azimuth(°),          P1(m),          P2(m),     L1(cycles),     L2(cycles),      Trop Delay,       Trop Map,  Relativity(m),      Sagnac(m), Tide Effect(m), Antenna Height, Sat Antenna(m),OffsetL1(cycles),OffsetL2(cycles),Windup(cycles)

Satellite Number: 10,GPS time：   91860,ztd:  2.3031,gps\_clock: -39412.36872110,glonass\_clock:      0.00000000

G19,0, -14617696.4714,  -9920046.2119, -20084205.0689,   -128075.6919,         4.7705,       136.4604,  25555726.1250,  25555734.4770, 134296194.3310, 104646453.0794,        24.1377,        11.2385,         4.6098,        23.0218,         0.0013,         0.0000,         0.8148,        -0.0154,         0.0154,         0.5224

G17,0, -18324889.6516,  15405540.0261,  11477464.8222,    -15739.1304,        14.5574,        24.6113,  24181809.1250,  24181820.5860, 127076286.7300,  99020528.3134,         8.9992,         3.9467,         6.1166,        12.8488,        -0.0071,         0.0000,         0.7958,         0.1128,         0.1212,         0.0501

G28,0, -12630049.0291,  18313959.0049, -13782301.2765,     95551.3318,        80.2511,        90.8865,  19732709.4610,  19732713.3010, 103696131.1610,  80802246.7935,         2.3367,         1.0146,         8.5468,         4.4298,        -0.0515,         0.0000,         1.0443,         0.4568,         0.4842,         0.2606

G07,0, -20443205.4026,  -1363756.9863, -17007860.8324,     88293.7244,        22.9805,       118.6209,  23365921.6330,  23365927.4060, 122788757.3870,  95679569.4565,         5.8582,         2.5536,        -4.4645,        25.0022,        -0.0133,         0.0000,         0.8274,         0.2045,         0.2007,         0.5551

G05,0,   3356903.5513,  24818435.1872,  -8644160.3943,   -118312.4341,        47.9401,       284.9715,  21493291.7270,  21493296.0940, 112947997.9020,  88011458.5495,         3.0988,         1.3463,         1.6878,       -18.2498,        -0.0469,         0.0000,         0.8207,         0.3880,         0.3872,        -0.7600
### 5.2成果数据文件格式
表头+解算的历元数+该历元解算出的X坐标+该历元解算出的Y坐标+该历元解算出的Z坐标+该历元解算出的X坐标改正数+该历元解算出的Y坐标改正数+该历元解算出的Z坐标改正数

例如：

1	-2364376.175498	4870295.955732	-3360814.971969	-38.790498	10.348832	-5.244069	

2	-2364334.848154	4870274.001229	-3360816.507765	2.536846	-11.605671	-6.779865	

3	-2364333.760433	4870275.193329	-3360808.371461	3.624567	-10.413571	1.356439	

4	-2364333.949337	4870278.252247	-3360805.321119	3.435663	-7.354653	4.406781	

5	-2364334.379239	4870281.054881	-3360804.642100	3.005761	-4.552019	5.085800	

6	-2364335.141636	4870284.090281	-3360805.490578	2.243364	-1.516619	4.237322	

7	-2364336.134308	4870286.302326	-3360806.725754	1.250692	0.695426	3.002146	

8	-2364337.036305	4870287.464707	-3360808.329688	0.348695	1.857807	1.398212	

9	-2364337.786921	4870287.965327	-3360809.865996	-0.401921	2.358427	-0.138096	

10	-2364338.418846	4870288.480393	-3360811.020078	-1.033846	2.873493	-1.292178

