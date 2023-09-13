#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QChartView>
#include <QLineSeries>
#include <QValueAxis>
#include <iostream>
#include<fstream>  //ifstream
#include <string>
#include <sstream>
#include <vector>
#include <typeinfo>
#include <math.h>
#include <Eigen/Dense>
#define PI acos(-1)

using namespace std;
using namespace Eigen;
QT_CHARTS_USE_NAMESPACE
typedef struct Data
{
    std::string PRN;
    double satnum, sats, satPos[3], sat_Clock, Elevation, Azimuth, Prange[2], L[2], Trop_Delay, Trop_Map, Relativity, Sagnac, Tide_Effect, Antenna_Height, Sat_Antenna, Tide_Offset[2], Windup;	// 城市坐标
}Data;//小结构体，用于存放单颗卫星的数据

typedef struct Satellite {
    int Satellite_Num;
    double gpst, ztd, gps_clock, glonass_clock;
    std::vector<Data> data;

}Satellite;//大结构体，用于存放整个历元的数据

double APPROX_POSITION[3];//全局变量
double X0=0, Y0=0, Z0=0;//全局变量存储接收机位置
MatrixXd POS(2, 2);//创建坐标结果矩阵
MatrixXd ANS(4, 1);
int geshu;

void txt_reader(std::string filename, std::vector<Satellite>& vec) {
    std::ifstream file;
    file.open(filename);
    std::string line;
    std::string temp_s;
    std::string temp[25];
    for (int i = 0; i < 11; i++) {//读取头文件
        std::getline(file, line);//一次一行，读取无用数据
        if (i == 4) {
            std::stringstream ss(line);
            int j = 0;
            while (std::getline(ss, temp_s, ',')) {
                if (temp_s.substr(temp_s.find_last_of(" "), line.size()).length() > 5)
                    APPROX_POSITION[j] = stod(temp_s.substr(temp_s.find_last_of(" "), line.size()));
                else
                    APPROX_POSITION[j] = stod(temp_s.substr(temp_s.find_first_of(" "), temp_s.find_last_of(" ")));
                j++;
            }
        }

    }
    while (std::getline(file, line)) {
        /*std::getline(file, line);*/
        Satellite satellite;
        std::stringstream ss(line);

        int i = 0;
        while (std::getline(ss, temp_s, ',')) {
            /*	std::cout << temp_s<<" "<< temp_s.find_last_of(" ") << " " << temp_s.find_last_of(":") <<" " << line.size() << std::endl;*/
            if (temp_s.find_last_of(" ") < line.size())
                temp[i] = temp_s.substr(temp_s.find_last_of(" "), line.size());
            else
                temp[i] = temp_s.substr(temp_s.find_last_of(":") + 1, line.size());

            i++;
        }
        ss.clear();
        satellite.Satellite_Num = stoi(temp[0]);
        satellite.gpst = stod(temp[1]);//string to double
        satellite.ztd = stod(temp[2]);
        satellite.gps_clock = stod(temp[3]);
        satellite.glonass_clock = stod(temp[4]);
        int j = 0;
        for (int i = 0; i < satellite.Satellite_Num; i++) {
            std::getline(file, line);
            std::stringstream ss(line);
            int j = 0;
            while (std::getline(ss, temp_s, ',')) {
                temp[j] = temp_s;
                j++;
            }

            Data data_temp;
            data_temp.PRN = temp[0];
            /*data_temp.satnum = stod(temp[1]);
            data_temp.sats = stod(temp[2]);*/
            data_temp.satPos[0] = stod(temp[2]);
            data_temp.satPos[1] = stod(temp[3]);
            data_temp.satPos[2] = stod(temp[4]);
            data_temp.sat_Clock = stod(temp[5]);
            data_temp.Elevation = stod(temp[6]);
            data_temp.Azimuth = stod(temp[7]);
            data_temp.Prange[0] = stod(temp[8]);
            data_temp.Prange[1] = stod(temp[9]);
            data_temp.L[0] = stod(temp[10]);
            data_temp.L[1] = stod(temp[11]);
            data_temp.Trop_Delay = stod(temp[12]);
            data_temp.Trop_Map = stod(temp[13]);
            data_temp.Relativity = stod(temp[14]);
            data_temp.Sagnac = stod(temp[15]);
            data_temp.Tide_Effect = stod(temp[16]);
            data_temp.Antenna_Height = stod(temp[17]);
            data_temp.Sat_Antenna = stod(temp[18]);
            data_temp.Tide_Offset[0] = stod(temp[19]);
            data_temp.Tide_Offset[1] = stod(temp[20]);
            data_temp.Windup = stod(temp[21]);
            satellite.data.push_back(data_temp);
            ss.clear();
        }
        vec.push_back(satellite);

    }
    //file.close();
}
std::vector<Satellite> satellite;//创建数据结构体


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->setWindowTitle("最优估计课程设计程序");
    setWindowIcon(QIcon(":/new/diqiu.ico"));
    //QString imgfile=QCoreApplication::applicationFilePath()+"lingmeng.jpg";
    //ui->showimage->setPixmap(QPixmap::fromImage(imgfile));
    //connect(ui->menu,&QMenu::triggered,this,&MainWindow::lineclear);
    QSound::play(":/new/cszx.wav");
    connect(ui->guanbi,&QAction::triggered,this,&MainWindow::guanbi);
    connect(ui->daoru,&QAction::triggered,this,&MainWindow::daoru);
    connect(ui->baocun,&QAction::triggered,this,&MainWindow::baocun);
    connect(ui->action_2,&QAction::triggered,this,&MainWindow::xinxi);
    connect(ui->Draw_X,&QAction::triggered,this,&MainWindow::Draw_X);
    connect(ui->Draw_Y,&QAction::triggered,this,&MainWindow::Draw_Y);
    connect(ui->Draw_Z,&QAction::triggered,this,&MainWindow::Draw_Z);

}

MainWindow::~MainWindow()
{
    //QString num3;

    delete ui;
}

void MainWindow::on_LSM_Button_clicked()
{
    ui->textBrowser->clear();
   //Qstring num2;
   //num1=ui->lineEdit->text();
   //double a=num1.toDouble();
   //cout<<a;
   //num2=QString::number(a,'f',2);
   //ui->textBrowser->setText(num2);
    MatrixXd zuobiao(3, 1);//将初始坐标存储到矩阵中
        for (int i = 0; i <= 2; i++)
        {
            zuobiao(i,0)= APPROX_POSITION[i];
            //cout << zuobiao(i, 0) << endl;
        }
        X0 = zuobiao(0, 0);//将初始坐标赋值给全局变量
        Y0 = zuobiao(1, 0);
        Z0 = zuobiao(2, 0);
        cout << endl;
        //MatrixXd H(1, 4);
        //double H1[100][4];
        //MatrixXd H(1, 4);
        MatrixXd H0(1, 4);//创建系数矩阵
        MatrixXd L(2, 2);//创建观测值矩阵
        //MatrixXd P(2, 2);
        //MatrixXd POS(2, 2);//创建坐标结果矩阵
        //MatrixXd LSM_ANS(4, 1);//创建改正数矩阵
        MatrixXd lsm_ans(4, 1);//创建改正数计算矩阵
        //MatrixXd ANS(2, 2);//创建结果矩阵
        MatrixXd P0(2, 2);//创建权矩阵
        MatrixXd fuzhu(2, 2);//创建辅助计算矩阵
        double dx;
        double dy;
        double dz;
        double R0;
        double l;
        double m;
        double n;
        //int geshu;
        int i = satellite.size();
        //qDebug()<<i;
        if(i==0)
        {
            QMessageBox::warning(this,"警告","请先导入文件");
        }
        QString geshu1;
        geshu1=ui->lineEdit->text();
        if(geshu1.isEmpty())
        {
            QMessageBox::warning(this,"警告","请输入需要解算的历元个数：");
        }
        else
        {
            geshu=geshu1.toInt();
            if (geshu>i||geshu<=0)//超出历元数据个数则结束程序
                {
                    //std::cout << "错误，超出最多解算个数" << endl;
                    QMessageBox::warning(this,"警告","错误，超出最多解算个数");
                    ui->textBrowser->insertPlainText("最多个数为：");
                    ui->textBrowser->insertPlainText(QString::number(i,'f',0));
                    ui->textBrowser->moveCursor(QTextCursor::End);
                    ui->textBrowser->append(QString(""));
                    //std::cout << "最多个数为：" << i << endl;
                    //return 0;
                }
                else
                {
                    //std::cout << "开始解算" << endl;
                    ui->textBrowser->insertPlainText("开始解算");
                    ui->textBrowser->moveCursor(QTextCursor::End);
                    ui->textBrowser->append(QString(""));
                    int Sat_num=0;
                        double Gps_clock = 0;
                        //MatrixXd::Constant (2, 4, 0.0);
                        //MatrixXd H(2, 4);
                        i = geshu;
                        POS.resize(3,geshu);//初始化定位结果矩阵
                        ANS.resize(4, geshu);//初始化改正数矩阵
                        for (int k = 0; k <= i - 1; k++)
                        {
                            Sat_num = satellite[k].Satellite_Num;//确定单个历元内的卫星个数
                            if (Sat_num <= 4)//若历元内的卫星数量小于4，则将改正数赋值为0，该次定位结果等于近似坐标，并跳过该历元
                            {
                                for (int a1 = 0; a1 <= 3; a1++)
                                {
                                    ANS(a1, k) = 0;
                                }
                                POS(0, k) = zuobiao(0);
                                POS(1, k) = zuobiao(1);
                                POS(2, k) = zuobiao(2);
                                continue;
                            }
                            else
                            {
                                MatrixXd H(Sat_num, 4);//初始化权阵
                                //MatrixXd P(2, 2);
                                //H.resize(Sat_num, 4);
                                L.resize(Sat_num, 1);//初始化观测值矩阵
                                //P.resize(Sat_num, Sat_num);
                                lsm_ans.resize(4, 1);//初始化改正数矩阵
                                //LSM_ANS.resize(4,1);
                                MatrixXd P = MatrixXd::Zero(Sat_num, Sat_num);//创建全0的以卫星数为行列数的方阵
                                /*for (int a1 = 0; a1 <= Sat_num-1; a1++)
                                {
                                    for (int a2 = 0; a2 <= Sat_num-1; a2++)
                                    {
                                        P(a1, a2) = 0;
                                    }
                                }*/
                                //P = MatrixXd::Zero();
                                Gps_clock = satellite[k].gps_clock;//提取该历元内的gps钟
                                for (int j = 0; j <= Sat_num - 1; j++)
                                {
                                    dx = satellite[k].data[j].satPos[0] - X0;//进行泰勒展开算法
                                    dy = satellite[k].data[j].satPos[1] - Y0;
                                    dz = satellite[k].data[j].satPos[2] - Z0;
                                    R0 = sqrt(pow(dx,2)+ pow(dy, 2) + pow(dz, 2));
                                    l = dx / R0;
                                    m = dy / R0;
                                    n = dz / R0;
                                    H0 << l, m, n, -1;
                                    //H.resize(j + 1, 4);
                                    H(j, 0) = l;//计算系数阵
                                    H(j, 1) = m;
                                    H(j, 2) = n;
                                    H(j, 3) = -1;
                                    L(j, 0) = satellite[k].data[j].Prange[0]-R0+satellite[k].data[j].sat_Clock-satellite[k].data[j].Trop_Delay-satellite[k].gps_clock;//计算改正数矩阵
                                    P(j, j) = pow(0.151, 2) / pow(sin(satellite[k].data[j].Elevation*PI/180), 2);//计算权矩阵
                                    //H(j, 0) = 1;
                                    //H(j, 1) = 2;
                                    //H(j, 2) = 3;
                                    //H(j, 3) = 4;
                                    //H << H,
                                    //	H0;
                                }
                                //cout << H << endl << endl;
                                 ui->textBrowser->insertPlainText("正在解算第");
                                 ui->textBrowser->insertPlainText(QString::number(k+1,'f',0));
                                 ui->textBrowser->insertPlainText("个历元");
                                 ui->textBrowser->moveCursor(QTextCursor::End);
                                 ui->textBrowser->append(QString(""));
                                 P = P.inverse();//计算逆矩阵
                                             H = -1 * H;
                                             fuzhu = H.transpose() * P * H;
                                             lsm_ans = fuzhu.inverse() * H.transpose() * P * L; //最小二乘法计算改正数
                                             X0 += lsm_ans(0);
                                             Y0 += lsm_ans(1);
                                             Z0 += lsm_ans(2);
                                             for (int a1 = 0; a1 <= 3; a1++)
                                             {
                                                 ANS(a1, k) = lsm_ans(a1);
                                                 ui->textBrowser->insertPlainText(QString::number(lsm_ans(a1),'f',6));
                                                 ui->textBrowser->insertPlainText("   ");
                                             }
                                             ui->textBrowser->moveCursor(QTextCursor::End);
                                             ui->textBrowser->append(QString(""));
                                             for (int a1 = 0; a1 <= 2; a1++)
                                             {
                                                  POS(a1, k) = zuobiao(a1, 0) + lsm_ans(a1);
                                                  //cout << POS(a1, k) << "   ";
                                                  ui->textBrowser->insertPlainText(QString::number(POS(a1, k),'f',6));
                                                  ui->textBrowser->insertPlainText("   ");
                                             }
                                             ui->textBrowser->moveCursor(QTextCursor::End);
                                             ui->textBrowser->append(QString(""));
                                             ui->textBrowser->moveCursor(QTextCursor::End);
                                             ui->textBrowser->append(QString(""));
                            }
                                }
                }
        }

        /*do
        {
            QMessageBox::warning(this,"警告","请输入需要解算的历元个数：");
            //geshu1=ui->lineEdit->text();

        }while(geshu1.isEmpty());
        */

}

void MainWindow::lineclear()
{
    ui->lineEdit->clear();
}


void MainWindow::guanbi()
{
    this->close();
}

void MainWindow::daoru()
{
    QString fileName = QFileDialog::getOpenFileName(this,"请选择一个文件",QCoreApplication::applicationFilePath(),"*.txt");
    if(fileName.isEmpty())
    {
        QMessageBox::warning(this,"警告","请选择一个文件");
    }
    else
    {
        qDebug()<<fileName;
    }
    string filename=fileName.toStdString();
    txt_reader(filename, satellite);
    ui->textBrowser->insertPlainText("接收机估计坐标位置：");
    ui->textBrowser->moveCursor(QTextCursor::End);
    ui->textBrowser->append(QString(""));
    ui->textBrowser->insertPlainText(QString::number(APPROX_POSITION[0],'f',4));
    ui->textBrowser->moveCursor(QTextCursor::End);
    ui->textBrowser->append(QString(""));
    ui->textBrowser->insertPlainText(QString::number(APPROX_POSITION[1],'f',4));
    ui->textBrowser->moveCursor(QTextCursor::End);
    ui->textBrowser->append(QString(""));
    ui->textBrowser->insertPlainText(QString::number(APPROX_POSITION[2],'f',4));
    ui->textBrowser->moveCursor(QTextCursor::End);
    ui->textBrowser->append(QString(""));
}

void MainWindow::xinxi()
{
    ui->textBrowser->clear();
    ui->textBrowser->insertPlainText("安徽理工大学");
    ui->textBrowser->moveCursor(QTextCursor::End);
    ui->textBrowser->append(QString(""));
    ui->textBrowser->insertPlainText("导航工程20-2");
    ui->textBrowser->moveCursor(QTextCursor::End);
    ui->textBrowser->append(QString(""));
    ui->textBrowser->insertPlainText("王杰");
    ui->textBrowser->moveCursor(QTextCursor::End);
    ui->textBrowser->append(QString(""));

}

void MainWindow::baocun()
{
    QString baocun=QFileDialog::getSaveFileName(this,"保存文件",QCoreApplication::applicationFilePath(),"*.txt");
    if(baocun.isEmpty())
    {
        QMessageBox::warning(this,"警告","请选择一个文件");
    }
    else
    {
        qDebug()<<baocun;
        string baocun2=baocun.toStdString();
        ofstream ans;
        ans.open(baocun2,ios::out);
        //ifstream file;
        //FILE *ans = file.open(baocun2);//fopen(baocun2, "w+");
        ans<<'\t'<<"X坐标"<<'\t'<<"Y坐标"<<'\t'<<"Z坐标"<<'\t'<<"X改正数"<<'\t'<<"Y改正数"<<'\t'<<"Z改正数"<<'\r'<<'\n';
        for (int a1=0;a1<=geshu-1;a1++)
        {
            //fprintf(ans, "%d\t", a1+1);
            ans<<a1+1<<'\t';
            for (int a2 = 0; a2 <= 2; a2++)
            {
                //fprintf(ans,"%f\t",POS(a2,a1));
                ans<<to_string(POS(a2,a1))<<'\t';

            }
            for (int a2 = 0; a2 <= 2; a2++)
            {
                //fprintf(ans, "%f\t", LSM_ANS(a2, a1));
                ans<<to_string(ANS(a2, a1))<<'\t';
            }
            //fprintf(ans,"\r\n");
            ans<<'\r'<<'\n';

        }
        ans.close();
        /*
            fprintf(ans, "\tX坐标\tY坐标\tZ坐标\tX改正数\tY改正数\tZ改正数\r\n");
            for (int a1=0;a1<=geshu-1;a1++)
            {
                fprintf(ans, "%d\t", a1+1);
                for (int a2 = 0; a2 <= 2; a2++)
                {
                    fprintf(ans,"%f\t",POS(a2,a1));
                }
                for (int a2 = 0; a2 <= 2; a2++)
                {
                    fprintf(ans, "%f\t", LSM_ANS(a2, a1));
                }
                fprintf(ans,"\r\n");

            }
            fclose(ans);
        */
    }

}

void MainWindow::on_KALMAN_Button_clicked()
{
    MatrixXd zuobiao(3, 1);//将初始坐标存储到矩阵中
        for (int i = 0; i <= 2; i++)
        {
            zuobiao(i, 0) = APPROX_POSITION[i];
            //cout << zuobiao(i, 0) << endl;
        }
        X0 = zuobiao(0, 0); Y0 = zuobiao(1, 0); Z0 = zuobiao(2, 0);
        //cout << endl;
        MatrixXd fai = MatrixXd::Zero(4, 4);
        for (int i = 0; i <= 3; i++)
        {
            fai(i, i) = 1;
        }
        MatrixXd Dx1 = MatrixXd::Zero(4, 4);
        Dx1(0, 0) = 25; Dx1(1, 1) = 25; Dx1(2, 2) = 25; Dx1(3, 3) = 10000;
        MatrixXd De = MatrixXd::Zero(4, 4);
        for (int i = 0; i <= 3; i++)
        {
            De(i, i) = 0.01;
        }
        MatrixXd H0(1, 4);
        MatrixXd H(1, 4);
        MatrixXd L(2, 2);
        //MatrixXd P(2, 2);
        //MatrixXd POS(2, 2);
        //MatrixXd LSM_ANS(4, 1);
        MatrixXd lsm_ans(4, 1);
        //MatrixXd ANS(2, 2);
        MatrixXd KALMAN_ANS(4, 1);
        MatrixXd kalman_ans(2, 2);
        MatrixXd X1(2, 2);
        MatrixXd X2(2, 2);
        MatrixXd Dx2(2, 2);
        MatrixXd Kk(2, 2);
        MatrixXd Vz(2, 2);
        MatrixXd fuzhu(2, 2);
        int k = 0;
        double dx = 0.0;
        double dy = 0.0;
        double dz = 0.0;
        double R0 = 0.0;
        double l = 0.0;
        double m = 0.0;
        double n = 0.0;
        int i = satellite.size();
        //int geshu;
        if(i==0)
        {
            QMessageBox::warning(this,"警告","请先导入文件");
        }
        QString geshu1;
        geshu1=ui->lineEdit->text();
         if(geshu1.isEmpty())
        {
            QMessageBox::warning(this,"警告","请输入需要解算的历元个数：");
        }
        //std::cout << "请输入需要解算的历元个数：" << endl;//输入解算历元的个数
        //std::cin >> geshu;
        else
        {
        //int i = satellite.size();
            geshu=geshu1.toInt();
            if (geshu > i || geshu <= 0)//超出历元数据个数则结束程序
            {
                //std::cout << "错误，超出最多解算个数" << endl;
                //std::cout << "最多个数为：" << i << endl;
                //return 0;
                QMessageBox::warning(this,"警告","错误，超出最多解算个数");
                ui->textBrowser->insertPlainText("最多个数为：");
                ui->textBrowser->insertPlainText(QString::number(i,'f',0));
                ui->textBrowser->moveCursor(QTextCursor::End);
                ui->textBrowser->append(QString(""));

            }
            else
            {
                //std::cout << "开始解算" << endl;
                ui->textBrowser->insertPlainText("开始解算");
                ui->textBrowser->moveCursor(QTextCursor::End);
                ui->textBrowser->append(QString(""));
                i = geshu;
                int a2 = 0;
                ANS.resize(4, i);
                POS.resize(3, i);
                int Sat_num = 0;
                double Gps_clock = 0;
                for (int k = 0; k <= i - 1; k++)
                {
                    Sat_num = satellite[k].Satellite_Num;
                    if (Sat_num <= 4)
                    {
                        for (int a1 = 0; a1 <= 3; a1++)
                        {
                            ANS(a1, k) = 0;
                        }
                        for (int a1 = 0; a1 <= 2; a1++)
                        {
                            POS(a1, k) = 0;
                        }
                        continue;
                    }
                    else
                    {
                        H.resize(Sat_num, 4);
                        L.resize(Sat_num, 1);
                        //P.resize(Sat_num, Sat_num);
                        lsm_ans.resize(4, 1);
                        MatrixXd P = MatrixXd::Zero(Sat_num, Sat_num);
                        Gps_clock = satellite[k].gps_clock;
                        for (int j = 0; j <= Sat_num - 1; j++)
                        {
                            dx = satellite[k].data[j].satPos[0] - X0;
                            dy = satellite[k].data[j].satPos[1] - Y0;
                            dz = satellite[k].data[j].satPos[2] - Z0;
                            R0 = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
                            l = dx / R0;
                            m = dy / R0;
                            n = dz / R0;
                            H0 << l, m, n, -1;
                            H(j, 0) = l;
                            H(j, 1) = m;
                            H(j, 2) = n;
                            H(j, 3) = -1;
                            L(j, 0) = satellite[k].data[j].Prange[0] - R0 + satellite[k].data[j].sat_Clock - satellite[k].data[j].Trop_Delay - satellite[k].gps_clock;
                            P(j, j) = pow(0.151, 2) / pow(sin(satellite[k].data[j].Elevation * PI / 180), 2);//要改
                        }
                        P = P.inverse();
                        H = -1 * H;
                        fuzhu = H.transpose() * P * H;
                        lsm_ans = fuzhu.inverse() * H.transpose() * P * L;
                    }
                    X0 += lsm_ans(0);
                    Y0 += lsm_ans(1);
                    Z0 += lsm_ans(2);
                    ANS(0, k) = lsm_ans(0);
                    ANS(1, k) = lsm_ans(1);
                    ANS(2, k) = lsm_ans(2);
                    ANS(3, k) = lsm_ans(3);
                    POS(0, k) = zuobiao(0) + lsm_ans(0);
                    POS(1, k) = zuobiao(1) + lsm_ans(1);
                    POS(2, k) = zuobiao(2) + lsm_ans(2);
                    a2 = k;
                    break;
                }
                kalman_ans = lsm_ans;
                for (int k = a2 + 1; k <= i - 1; k++)
                {
                    Sat_num = satellite[k].Satellite_Num;
                    if (Sat_num <= 4)
                    {
                        for (int a1 = 0; a1 <= 3; a1++)
                        {
                            ANS(a1, k) = 0;
                        }
                        POS(0, k) = zuobiao(0);
                        POS(1, k) = zuobiao(1);
                        POS(2, k) = zuobiao(2);
                        continue;
                    }
                    else
                    {
                        H.resize(Sat_num, 4);
                        L.resize(Sat_num, 1);
                        //kalman_ans.resize(4, 1);
                        Gps_clock = satellite[k].gps_clock;
                        MatrixXd P = MatrixXd::Zero(Sat_num, Sat_num);
                        for (int j = 0; j <= Sat_num - 1; j++)
                        {
                            dx = satellite[k].data[j].satPos[0] - X0;
                            dy = satellite[k].data[j].satPos[1] - Y0;
                            dz = satellite[k].data[j].satPos[2] - Z0;
                            R0 = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
                            l = dx / R0;
                            m = dy / R0;
                            n = dz / R0;
                            H0 << l, m, n, -1;
                            H(j, 0) = l;
                            H(j, 1) = m;
                            H(j, 2) = n;
                            H(j, 3) = -1;
                            L(j, 0) = satellite[k].data[j].Prange[0] - R0 + satellite[k].data[j].sat_Clock - satellite[k].data[j].Trop_Delay - satellite[k].gps_clock;
                            P(j, j) = pow(0.151, 2) / pow(sin(satellite[k].data[j].Elevation * PI / 180), 2);//要改
                        }
                        ui->textBrowser->insertPlainText("正在解算第");
                        ui->textBrowser->insertPlainText(QString::number(k+1,'f',0));
                        ui->textBrowser->insertPlainText("个历元");
                        ui->textBrowser->moveCursor(QTextCursor::End);
                        ui->textBrowser->append(QString(""));
                        //cout << "正在解算第" << k + 1 << "个历元" << endl;
                        P = P.inverse();
                        H = -1 * H;
                        MatrixXd X1(2, 2);
                        X1 = fai * kalman_ans;//要改
                        Dx2 = fai * Dx1 * fai.transpose() + De;
                        fuzhu = H * Dx2 * H.transpose() + P;
                        Kk = Dx2 * H.transpose() * fuzhu.inverse();
                        Vz = L - H * X1;
                        X2 = X1 + Kk * Vz;
                        kalman_ans = X2;
                        fuzhu = Kk * H;
                        int fuzhu1 = sqrt(fuzhu.size());
                        MatrixXd eyezhen = MatrixXd::Identity(fuzhu1,fuzhu1);
                        Dx1 = (eyezhen - fuzhu) * Dx2;
                        ANS(0, k) = kalman_ans(0);
                        ANS(1, k) = kalman_ans(1);
                        ANS(2, k) = kalman_ans(2);
                        ANS(3, k) = kalman_ans(3);
                        X0 += kalman_ans(0);
                        Y0 += kalman_ans(1);
                        Z0 += kalman_ans(2);
                        //POS(0, k) = zuobiao(0) + kalman_ans(0);
                        //POS(1, k) = zuobiao(1) + kalman_ans(1);
                        //POS(2, k) = zuobiao(2) + kalman_ans(2);
                        for (int a1 = 0; a1 <= 3; a1++)
                        {
                            //ANS(a1, k) = lsm_ans(a1);
                            ui->textBrowser->insertPlainText(QString::number(ANS(a1, k),'f',6));
                            ui->textBrowser->insertPlainText("   ");
                        }
                        ui->textBrowser->moveCursor(QTextCursor::End);
                        ui->textBrowser->append(QString(""));
                        for (int a1 = 0; a1 <= 2; a1++)
                        {
                             POS(a1, k) = zuobiao(a1, 0) + kalman_ans(a1);
                             //cout << POS(a1, k) << "   ";
                             ui->textBrowser->insertPlainText(QString::number(POS(a1, k),'f',6));
                             ui->textBrowser->insertPlainText("   ");
                        }
                        ui->textBrowser->moveCursor(QTextCursor::End);
                        ui->textBrowser->append(QString(""));
                        ui->textBrowser->moveCursor(QTextCursor::End);
                        ui->textBrowser->append(QString(""));
                        kalman_ans.resize(4, 1);
                        //cout << POS(0, k) << endl << POS(1, k) << endl << POS(2, k) << endl << endl;
                        //cout << KALMAN_ANS(0, k) << endl << KALMAN_ANS(1, k) << endl << KALMAN_ANS(2, k) << endl << KALMAN_ANS(3, k) << endl << endl;
                        //cout << eyezhen;
                        //cout << fuzhu << endl << endl;
                        //cout << fuzhu.size();
                        //fuzhu.size;

                    }

                 }
            }
        }


}

void MainWindow::on_KALMANCVButton_clicked()
{
    MatrixXd zuobiao(3, 1);//将初始坐标存储到矩阵中
        for (int i = 0; i <= 2; i++)
        {
            zuobiao(i, 0) = APPROX_POSITION[i];
            //cout << zuobiao(i, 0) << endl;
        }
        X0 = zuobiao(0, 0); Y0 = zuobiao(1, 0); Z0 = zuobiao(2, 0);
        //cout << endl;
        MatrixXd fai = MatrixXd::Zero(7, 7);
        for (int i = 0; i <= 6; i++)
        {
            fai(i, i) = 1;
        }
        for (int i = 0; i <= 2; i++)
        {
            fai(i, i+3) = 30;
        }
        MatrixXd Dx1 = MatrixXd::Zero(7, 7);
        Dx1(0, 0) = 25; Dx1(1, 1) = 25; Dx1(2, 2) = 25; Dx1(3, 3) = 1;Dx1(4, 4) = 1;Dx1(5, 5) = 1;Dx1(6, 6) = 10000;
        MatrixXd De = MatrixXd::Zero(7, 7);
        for (int i = 0; i <= 6; i++)
        {
            De(i, i) = 0.01;
        }
        MatrixXd H0(1, 4);
        MatrixXd H(1, 4);
        MatrixXd L(2, 2);
        //MatrixXd P(2, 2);
        //MatrixXd POS(2, 2);
        //MatrixXd LSM_ANS(4, 1);
        MatrixXd lsm_ans(4, 1);
        //MatrixXd ANS(2, 2);
        MatrixXd KALMAN_ANS(4, 1);
        MatrixXd kalman_ans(7, 1);
        MatrixXd X1(2, 2);
        MatrixXd X2(2, 2);
        MatrixXd Dx2(2, 2);
        MatrixXd Kk(2, 2);
        MatrixXd Vz(2, 2);
        MatrixXd fuzhu(2, 2);
        int k = 0;
        double dx = 0.0;
        double dy = 0.0;
        double dz = 0.0;
        double R0 = 0.0;
        double l = 0.0;
        double m = 0.0;
        double n = 0.0;
        int i = satellite.size();
        //int geshu;
        if(i==0)
        {
            QMessageBox::warning(this,"警告","请先导入文件");
        }
        QString geshu1;
        geshu1=ui->lineEdit->text();
         if(geshu1.isEmpty())
        {
            QMessageBox::warning(this,"警告","请输入需要解算的历元个数：");
        }
        //std::cout << "请输入需要解算的历元个数：" << endl;//输入解算历元的个数
        //std::cin >> geshu;
        else
        {
        //int i = satellite.size();
            geshu=geshu1.toInt();
            if (geshu > i || geshu <= 0)//超出历元数据个数则结束程序
            {
                //std::cout << "错误，超出最多解算个数" << endl;
                //std::cout << "最多个数为：" << i << endl;
                //return 0;
                QMessageBox::warning(this,"警告","错误，超出最多解算个数");
                ui->textBrowser->insertPlainText("最多个数为：");
                ui->textBrowser->insertPlainText(QString::number(i,'f',0));
                ui->textBrowser->moveCursor(QTextCursor::End);
                ui->textBrowser->append(QString(""));

            }
            else
            {
                //std::cout << "开始解算" << endl;
                ui->textBrowser->insertPlainText("开始解算");
                ui->textBrowser->moveCursor(QTextCursor::End);
                ui->textBrowser->append(QString(""));
                i = geshu;
                int a2 = 0;
                ANS.resize(4, i);
                POS.resize(3, i);
                int Sat_num = 0;
                double Gps_clock = 0;
                for (int k = 0; k <= i - 1; k++)
                {
                    Sat_num = satellite[k].Satellite_Num;
                    if (Sat_num <= 4)
                    {
                        for (int a1 = 0; a1 <= 3; a1++)
                        {
                            ANS(a1, k) = 0;
                        }
                        for (int a1 = 0; a1 <= 2; a1++)
                        {
                            POS(a1, k) = 0;
                        }
                        continue;
                    }
                    else
                    {
                        H.resize(Sat_num, 4);
                        L.resize(Sat_num, 1);
                        //P.resize(Sat_num, Sat_num);
                        lsm_ans.resize(4, 1);
                        MatrixXd P = MatrixXd::Zero(Sat_num, Sat_num);
                        Gps_clock = satellite[k].gps_clock;
                        for (int j = 0; j <= Sat_num - 1; j++)
                        {
                            dx = satellite[k].data[j].satPos[0] - X0;
                            dy = satellite[k].data[j].satPos[1] - Y0;
                            dz = satellite[k].data[j].satPos[2] - Z0;
                            R0 = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
                            l = dx / R0;
                            m = dy / R0;
                            n = dz / R0;
                            H0 << l, m, n, -1;
                            H(j, 0) = l;
                            H(j, 1) = m;
                            H(j, 2) = n;
                            H(j, 3) = -1;
                            L(j, 0) = satellite[k].data[j].Prange[0] - R0 + satellite[k].data[j].sat_Clock - satellite[k].data[j].Trop_Delay - satellite[k].gps_clock;
                            P(j, j) = pow(0.151, 2) / pow(sin(satellite[k].data[j].Elevation * PI / 180), 2);//要改
                        }
                        P = P.inverse();
                        H = -1 * H;
                        fuzhu = H.transpose() * P * H;
                        lsm_ans = fuzhu.inverse() * H.transpose() * P * L;
                    }
                    X0 += lsm_ans(0);
                    Y0 += lsm_ans(1);
                    Z0 += lsm_ans(2);
                    ANS(0, k) = lsm_ans(0);
                    ANS(1, k) = lsm_ans(1);
                    ANS(2, k) = lsm_ans(2);
                    ANS(3, k) = lsm_ans(3);
                    POS(0, k) = zuobiao(0) + lsm_ans(0);
                    POS(1, k) = zuobiao(1) + lsm_ans(1);
                    POS(2, k) = zuobiao(2) + lsm_ans(2);
                    a2 = k;
                    break;
                }
                //kalman_ans = lsm_ans;
                for(int a3=0;a3<=2;a3++)
                {
                    kalman_ans(a3)=lsm_ans(a3);
                }
                kalman_ans(6)=lsm_ans(3);
                for(int a3=3;a3<=5;a3++)
                {
                    kalman_ans(a3)=0;
                }
                for (int k = a2 + 1; k <= i - 1; k++)
                {
                    Sat_num = satellite[k].Satellite_Num;
                    if (Sat_num <= 4)
                    {
                        for (int a1 = 0; a1 <= 3; a1++)
                        {
                            ANS(a1, k) = 0;
                        }
                        POS(0, k) = zuobiao(0);
                        POS(1, k) = zuobiao(1);
                        POS(2, k) = zuobiao(2);
                        continue;
                    }
                    else
                    {
                        H.resize(Sat_num, 7);
                        L.resize(Sat_num, 1);
                        //kalman_ans.resize(4, 1);
                        Gps_clock = satellite[k].gps_clock;
                        MatrixXd P = MatrixXd::Zero(Sat_num, Sat_num);
                        for (int j = 0; j <= Sat_num - 1; j++)
                        {
                            dx = satellite[k].data[j].satPos[0] - X0;
                            dy = satellite[k].data[j].satPos[1] - Y0;
                            dz = satellite[k].data[j].satPos[2] - Z0;
                            R0 = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
                            l = dx / R0;
                            m = dy / R0;
                            n = dz / R0;
                            H0 << l, m, n, -1;
                            H(j, 0) = l;
                            H(j, 1) = m;
                            H(j, 2) = n;
                            H(j, 3) = 0;
                            H(j, 4) = 0;
                            H(j, 5) = 0;
                            H(j, 6) = -1;
                            L(j, 0) = satellite[k].data[j].Prange[0] - R0 + satellite[k].data[j].sat_Clock - satellite[k].data[j].Trop_Delay - satellite[k].gps_clock;
                            P(j, j) = pow(0.151, 2) / pow(sin(satellite[k].data[j].Elevation * PI / 180), 2);//要改
                        }
                        ui->textBrowser->insertPlainText("正在解算第");
                        ui->textBrowser->insertPlainText(QString::number(k+1,'f',0));
                        ui->textBrowser->insertPlainText("个历元");
                        ui->textBrowser->moveCursor(QTextCursor::End);
                        ui->textBrowser->append(QString(""));
                        //cout << "正在解算第" << k + 1 << "个历元" << endl;
                        P = P.inverse();
                        H = -1 * H;
                        MatrixXd X1(2, 2);
                        X1 = fai * kalman_ans;//要改
                        Dx2 = fai * Dx1 * fai.transpose() + De;
                        fuzhu = H * Dx2 * H.transpose() + P;
                        Kk = Dx2 * H.transpose() * fuzhu.inverse();
                        Vz = L - H * X1;
                        X2 = X1 + Kk * Vz;
                        kalman_ans = X2;
                        fuzhu = Kk * H;
                        int fuzhu1 = sqrt(fuzhu.size());
                        MatrixXd eyezhen = MatrixXd::Identity(fuzhu1,fuzhu1);
                        Dx1 = (eyezhen - fuzhu) * Dx2;
                        ANS(0, k) = kalman_ans(0);
                        ANS(1, k) = kalman_ans(1);
                        ANS(2, k) = kalman_ans(2);
                        ANS(3, k) = kalman_ans(6);
                        X0 += kalman_ans(0);
                        Y0 += kalman_ans(1);
                        Z0 += kalman_ans(2);
                        //POS(0, k) = zuobiao(0) + kalman_ans(0);
                        //POS(1, k) = zuobiao(1) + kalman_ans(1);
                        //POS(2, k) = zuobiao(2) + kalman_ans(2);
                        for (int a1 = 0; a1 <= 3; a1++)
                        {
                            //ANS(a1, k) = lsm_ans(a1);
                            ui->textBrowser->insertPlainText(QString::number(ANS(a1, k),'f',6));
                            ui->textBrowser->insertPlainText("   ");
                        }
                        ui->textBrowser->moveCursor(QTextCursor::End);
                        ui->textBrowser->append(QString(""));
                        for (int a1 = 0; a1 <= 2; a1++)
                        {
                             POS(a1, k) = zuobiao(a1, 0) + kalman_ans(a1);
                             //cout << POS(a1, k) << "   ";
                             ui->textBrowser->insertPlainText(QString::number(POS(a1, k),'f',6));
                             ui->textBrowser->insertPlainText("   ");
                        }
                        ui->textBrowser->moveCursor(QTextCursor::End);
                        ui->textBrowser->append(QString(""));
                        ui->textBrowser->moveCursor(QTextCursor::End);
                        ui->textBrowser->append(QString(""));
                        kalman_ans.resize( 7, 1);
                        //cout << POS(0, k) << endl << POS(1, k) << endl << POS(2, k) << endl << endl;
                        //cout << KALMAN_ANS(0, k) << endl << KALMAN_ANS(1, k) << endl << KALMAN_ANS(2, k) << endl << KALMAN_ANS(3, k) << endl << endl;
                        //cout << eyezhen;
                        //cout << fuzhu << endl << endl;
                        //cout << fuzhu.size();
                        //fuzhu.size;

                    }

                 }
            }
        }
}

void MainWindow::on_KALMANCAButton_clicked()
{
    MatrixXd zuobiao(3, 1);//将初始坐标存储到矩阵中
        for (int i = 0; i <= 2; i++)
        {
            zuobiao(i, 0) = APPROX_POSITION[i];
            //cout << zuobiao(i, 0) << endl;
        }
        X0 = zuobiao(0, 0); Y0 = zuobiao(1, 0); Z0 = zuobiao(2, 0);
        //cout << endl;
        MatrixXd fai = MatrixXd::Zero(10, 10);
        for (int i = 0; i <= 9; i++)
        {
            fai(i, i) = 1;
        }
        for (int i = 0; i <= 5; i++)
        {
            fai(i, i+3) = 30;
        }
        for (int i = 0; i <= 2; i++)
        {
            fai(i, i+6) = 450;
        }
        MatrixXd Dx1 = MatrixXd::Zero(10, 10);
        Dx1(0, 0) = 25; Dx1(1, 1) = 25; Dx1(2, 2) = 25; Dx1(3, 3) = 1;Dx1(4, 4) = 1;Dx1(5, 5) = 1;Dx1(6, 6) = 1;Dx1(7, 7) = 1;Dx1(8, 8) = 1;Dx1(9, 9) = 10000;
        MatrixXd De = MatrixXd::Zero(10, 10);
        for (int i = 0; i <=9; i++)
        {
            De(i, i) = 0.01;
        }
        MatrixXd H0(1, 4);
        MatrixXd H(1, 4);
        MatrixXd L(2, 2);
        //MatrixXd P(2, 2);
        //MatrixXd POS(2, 2);
        //MatrixXd LSM_ANS(4, 1);
        MatrixXd lsm_ans(4, 1);
        //MatrixXd ANS(2, 2);
        MatrixXd KALMAN_ANS(4, 1);
        MatrixXd kalman_ans(10, 1);
        MatrixXd X1(2, 2);
        MatrixXd X2(2, 2);
        MatrixXd Dx2(2, 2);
        MatrixXd Kk(2, 2);
        MatrixXd Vz(2, 2);
        MatrixXd fuzhu(2, 2);
        int k = 0;
        double dx = 0.0;
        double dy = 0.0;
        double dz = 0.0;
        double R0 = 0.0;
        double l = 0.0;
        double m = 0.0;
        double n = 0.0;
        int i = satellite.size();
        //int geshu;
        if(i==0)
        {
            QMessageBox::warning(this,"警告","请先导入文件");
        }
        QString geshu1;
        geshu1=ui->lineEdit->text();
         if(geshu1.isEmpty())
        {
            QMessageBox::warning(this,"警告","请输入需要解算的历元个数：");
        }
        //std::cout << "请输入需要解算的历元个数：" << endl;//输入解算历元的个数
        //std::cin >> geshu;
        else
        {
        //int i = satellite.size();
            geshu=geshu1.toInt();
            if (geshu > i || geshu <= 0)//超出历元数据个数则结束程序
            {
                //std::cout << "错误，超出最多解算个数" << endl;
                //std::cout << "最多个数为：" << i << endl;
                //return 0;
                QMessageBox::warning(this,"警告","错误，超出最多解算个数");
                ui->textBrowser->insertPlainText("最多个数为：");
                ui->textBrowser->insertPlainText(QString::number(i,'f',0));
                ui->textBrowser->moveCursor(QTextCursor::End);
                ui->textBrowser->append(QString(""));

            }
            else
            {
                //std::cout << "开始解算" << endl;
                ui->textBrowser->insertPlainText("开始解算");
                ui->textBrowser->moveCursor(QTextCursor::End);
                ui->textBrowser->append(QString(""));
                i = geshu;
                int a2 = 0;
                ANS.resize(4, i);
                POS.resize(3, i);
                int Sat_num = 0;
                double Gps_clock = 0;
                for (int k = 0; k <= i - 1; k++)
                {
                    Sat_num = satellite[k].Satellite_Num;
                    if (Sat_num <= 4)
                    {
                        for (int a1 = 0; a1 <= 3; a1++)
                        {
                            ANS(a1, k) = 0;
                        }
                        for (int a1 = 0; a1 <= 2; a1++)
                        {
                            POS(a1, k) = 0;
                        }
                        continue;
                    }
                    else
                    {
                        H.resize(Sat_num, 4);
                        L.resize(Sat_num, 1);
                        //P.resize(Sat_num, Sat_num);
                        lsm_ans.resize(4, 1);
                        MatrixXd P = MatrixXd::Zero(Sat_num, Sat_num);
                        Gps_clock = satellite[k].gps_clock;
                        for (int j = 0; j <= Sat_num - 1; j++)
                        {
                            dx = satellite[k].data[j].satPos[0] - X0;
                            dy = satellite[k].data[j].satPos[1] - Y0;
                            dz = satellite[k].data[j].satPos[2] - Z0;
                            R0 = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
                            l = dx / R0;
                            m = dy / R0;
                            n = dz / R0;
                            H0 << l, m, n, -1;
                            H(j, 0) = l;
                            H(j, 1) = m;
                            H(j, 2) = n;
                            H(j, 3) = -1;
                            L(j, 0) = satellite[k].data[j].Prange[0] - R0 + satellite[k].data[j].sat_Clock - satellite[k].data[j].Trop_Delay - satellite[k].gps_clock;
                            P(j, j) = pow(0.151, 2) / pow(sin(satellite[k].data[j].Elevation * PI / 180), 2);//要改
                        }
                        P = P.inverse();
                        H = -1 * H;
                        fuzhu = H.transpose() * P * H;
                        lsm_ans = fuzhu.inverse() * H.transpose() * P * L;
                    }
                    X0 += lsm_ans(0);
                    Y0 += lsm_ans(1);
                    Z0 += lsm_ans(2);
                    ANS(0, k) = lsm_ans(0);
                    ANS(1, k) = lsm_ans(1);
                    ANS(2, k) = lsm_ans(2);
                    ANS(3, k) = lsm_ans(3);
                    POS(0, k) = zuobiao(0) + lsm_ans(0);
                    POS(1, k) = zuobiao(1) + lsm_ans(1);
                    POS(2, k) = zuobiao(2) + lsm_ans(2);
                    a2 = k;
                    break;
                }
                //kalman_ans = lsm_ans;
                for(int a3=0;a3<=2;a3++)
                {
                    kalman_ans(a3)=lsm_ans(a3);
                }
                kalman_ans(9)=lsm_ans(3);
                for(int a3=3;a3<=8;a3++)
                {
                    kalman_ans(a3)=0;
                }
                for (int k = a2 + 1; k <= i - 1; k++)
                {
                    Sat_num = satellite[k].Satellite_Num;
                    if (Sat_num <= 4)
                    {
                        for (int a1 = 0; a1 <= 3; a1++)
                        {
                            ANS(a1, k) = 0;
                        }
                        POS(0, k) = zuobiao(0);
                        POS(1, k) = zuobiao(1);
                        POS(2, k) = zuobiao(2);
                        continue;
                    }
                    else
                    {
                        H.resize(Sat_num, 10);
                        L.resize(Sat_num, 1);
                        //kalman_ans.resize(4, 1);
                        Gps_clock = satellite[k].gps_clock;
                        MatrixXd P = MatrixXd::Zero(Sat_num, Sat_num);
                        for (int j = 0; j <= Sat_num - 1; j++)
                        {
                            dx = satellite[k].data[j].satPos[0] - X0;
                            dy = satellite[k].data[j].satPos[1] - Y0;
                            dz = satellite[k].data[j].satPos[2] - Z0;
                            R0 = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
                            l = dx / R0;
                            m = dy / R0;
                            n = dz / R0;
                            H0 << l, m, n, -1;
                            H(j, 0) = l;
                            H(j, 1) = m;
                            H(j, 2) = n;
                            H(j, 3) = 0;
                            H(j, 4) = 0;
                            H(j, 5) = 0;
                            H(j, 6) = 0;
                            H(j, 7) = 0;
                            H(j, 8) = 0;
                            H(j, 9) = -1;
                            L(j, 0) = satellite[k].data[j].Prange[0] - R0 + satellite[k].data[j].sat_Clock - satellite[k].data[j].Trop_Delay - satellite[k].gps_clock;
                            P(j, j) = pow(0.151, 2) / pow(sin(satellite[k].data[j].Elevation * PI / 180), 2);//要改
                        }
                        ui->textBrowser->insertPlainText("正在解算第");
                        ui->textBrowser->insertPlainText(QString::number(k+1,'f',0));
                        ui->textBrowser->insertPlainText("个历元");
                        ui->textBrowser->moveCursor(QTextCursor::End);
                        ui->textBrowser->append(QString(""));
                        //cout << "正在解算第" << k + 1 << "个历元" << endl;
                        P = P.inverse();
                        H = -1 * H;
                        MatrixXd X1(2, 2);
                        X1 = fai * kalman_ans;//要改
                        Dx2 = fai * Dx1 * fai.transpose() + De;
                        fuzhu = H * Dx2 * H.transpose() + P;
                        Kk = Dx2 * H.transpose() * fuzhu.inverse();
                        Vz = L - H * X1;
                        X2 = X1 + Kk * Vz;
                        kalman_ans = X2;
                        fuzhu = Kk * H;
                        int fuzhu1 = sqrt(fuzhu.size());
                        MatrixXd eyezhen = MatrixXd::Identity(fuzhu1,fuzhu1);
                        Dx1 = (eyezhen - fuzhu) * Dx2;
                        ANS(0, k) = kalman_ans(0);
                        ANS(1, k) = kalman_ans(1);
                        ANS(2, k) = kalman_ans(2);
                        ANS(3, k) = kalman_ans(9);
                        X0 += kalman_ans(0);
                        Y0 += kalman_ans(1);
                        Z0 += kalman_ans(2);
                        //POS(0, k) = zuobiao(0) + kalman_ans(0);
                        //POS(1, k) = zuobiao(1) + kalman_ans(1);
                        //POS(2, k) = zuobiao(2) + kalman_ans(2);
                        for (int a1 = 0; a1 <= 3; a1++)
                        {
                            //ANS(a1, k) = lsm_ans(a1);
                            ui->textBrowser->insertPlainText(QString::number(ANS(a1, k),'f',6));
                            ui->textBrowser->insertPlainText("   ");
                        }
                        ui->textBrowser->moveCursor(QTextCursor::End);
                        ui->textBrowser->append(QString(""));
                        for (int a1 = 0; a1 <= 2; a1++)
                        {
                             POS(a1, k) = zuobiao(a1, 0) + kalman_ans(a1);
                             //cout << POS(a1, k) << "   ";
                             ui->textBrowser->insertPlainText(QString::number(POS(a1, k),'f',6));
                             ui->textBrowser->insertPlainText("   ");
                        }
                        ui->textBrowser->moveCursor(QTextCursor::End);
                        ui->textBrowser->append(QString(""));
                        ui->textBrowser->moveCursor(QTextCursor::End);
                        ui->textBrowser->append(QString(""));
                        kalman_ans.resize(10, 1);
                        //cout << POS(0, k) << endl << POS(1, k) << endl << POS(2, k) << endl << endl;
                        //cout << KALMAN_ANS(0, k) << endl << KALMAN_ANS(1, k) << endl << KALMAN_ANS(2, k) << endl << KALMAN_ANS(3, k) << endl << endl;
                        //cout << eyezhen;
                        //cout << fuzhu << endl << endl;
                        //cout << fuzhu.size();
                        //fuzhu.size;

                    }

                 }
            }
        }
}

void MainWindow::Draw_X()
{
    if(POS.size()==4)
    {
        QMessageBox::warning(this,"警告","请先进行解算");
    }
    else
    {
    drawxx *drawxxw=new drawxx(POS,geshu,APPROX_POSITION);
    //Form *formw = new Form(POS,geshu,APPROX_POSITION);
    //formw=new Form;
    drawxxw->show();
    }
    /*
     *
    QChartView* xdraw=new QChartView(this);
    QChart* chart=new QChart();
    xdraw->setChart(chart);
    setCentralWidget(xdraw);
    //formw->show();

    QLineSeries* Drawx=new QLineSeries;
    Drawx->setName("X坐标");
    chart->addSeries(Drawx);
    for(int a1=0;a1<=geshu-1;a1++)
    {
        Drawx->append(a1,POS(0,a1));
    }
    QValueAxis* axisX=new QValueAxis;
    axisX->setRange(1,geshu);
    chart->setAxisX(axisX,Drawx);

    QValueAxis* axisY=new QValueAxis;
    axisY->setRange(APPROX_POSITION[0]+15,APPROX_POSITION[0]-15);
    chart->setAxisY(axisY,Drawx);
*/

    //QChartView* xdraw=new QChartView(this);
    //QChart* chart=new QChart();


}

void MainWindow::Draw_Y()
{

    if(POS.size()==4)
    {
        QMessageBox::warning(this,"警告","请先进行解算");
    }
    else
    {
    drawyy *drawyyw=new drawyy(POS,geshu,APPROX_POSITION);
    //Form *formw = new Form(POS,geshu,APPROX_POSITION);
    //formw=new Form;
    drawyyw->show();
    }
}

void MainWindow::Draw_Z()
{
    if(POS.size()==4)
    {
        QMessageBox::warning(this,"警告","请先进行解算");
    }
    else
    {
    drawzz *drawzzw=new drawzz(POS,geshu,APPROX_POSITION);
    //Form *formw = new Form(POS,geshu,APPROX_POSITION);
    //formw=new Form;
    drawzzw->show();
    }
}
