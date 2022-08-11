#ifndef PLandEE_H
#define PLandEE_H
#include<vector>
#include<cstdio>
#include<iostream>
#include<cmath>
#include<ctime>
#include<algorithm>
#include "Point.h"
#include "Kmeans.h"

using namespace std;
//数据转换
vector<double> GeographicToGeodeticRectangular(double Long, double Lati, double Height);
vector<double> GeodeticRectangularToStationRectangular(double xg, double yg, double zg, double Long, double Lati, double Height);

vector<double> StationRectangularToGeodeticRectangular(double xStationTarget, double yStationTarget, double zStationTarget, double longitudeStation, double latitudeStation, double heightStation);
vector<double> GeodeticRectangularToGeographic(double dX, double dY, double dZ);

double PicAngle2AlgAngle(double in);
double AlgAngle2PicAngle(double in);

// DF data 测向数据 也是最原始的数据
class DFdata{
private:
	double x; //经度
	double y; //纬度
	double z; //高度
	double alpha; //平面的偏角  0  --- 2pi
	double beta;  //高度偏角   -pi --- pi
	double r[2];	  // 当前坐标与所设原点的距离
	double sk[2];	  // 测量的从测量点到定位点的单位向量
	double distance; // 与定位点的距离 xy平面上的
public:
	DFdata();
	DFdata(double x_in, double y_in, double z_in, double alpha_in, double beta_in, bool isPI); // in  input的意思  

	//属性 获取和更改
	void set_x(double x_in);
	void set_y(double y_in);
	void set_z(double z_in);
	double set_distance(double x_in, double y_in);	//算出的定位点 到 测量点的距离
	double get_distance();
	// 关于角度这里的设置， 角度有两个单位 弧度 和 角度 传入的时候说明是弧度还是角度 isPI是弧度  isPI = 1 弧度  isPI = 0 角度
	void set_alpha(double alpha_in, bool isPI);
	void set_beta(double beta_in, bool isPI);
	double get_alpha();
	double get_beta();
	vector<double> get_point();  //3*1
	// 有关算法中参数r的式子  r 2*1 二维向量
	vector<double> set_r(double x_in, double y_in);
	vector<double> get_r();
	// 向量sk相关函数 sk 2*1 二维向量 sk原本是从测量点 到 定位点垂直于测向线上的点， 但是因为定位点不固定，我们这里用测向线方向的单位向量
	vector<double> get_sk();
	vector<double> get_skt();
	// 判断是否反向  射线 非直线
	bool isValid(double x, double y);
};

// 基于测量数据的无源定位解法  Q1 算法一
class PLocation{
private:
	vector<DFdata> dflist; // 同一批次的原始数据
	vector<vector<double>> A;		// 矩阵A skt  N*2 
	vector<double> b;		// 矩阵b   skt * r  N*1
	int n;					//这一批次的数据的数量
	double x;					//定位点经度
	double y;					//定位点维度
	double z;					//定位点高度
	double Pe;					//误差椭圆准确度
	double rl;					//长轴 r long
	double rs;					//短轴 r short
	double fan;					//偏角
	double sig11;
	double sig12;
	double sig22;
	bool Valid;					//是否有效
	vector<double> stationloc;	//辅助站心
	vector<vector<double>> result;		//输出结果
	//辅助当前算法函数值
	vector<vector<double>> pinvA; //pinv(A) 矩阵 (ATA)^-1 * A    2 * N
	double help_b(vector<double> skt, vector<double> r); //辅助生成b的参数
	// 算法的执行步骤
	vector<vector<double>> get_ATA(); //获得AT * A  2*2 的矩阵
	bool ATAsuit(vector<vector<double>> ATA); //判断ATA是否合法
	vector<vector<double>> get_reATA(vector<vector<double>> ATA);	//	获取ATA的逆
	vector<vector<double>> get_pinvA();//获取算法中参数pinv(A)  即 (ATA)^-1 * A    2 * N的矩阵
	vector<double> set_locxy(); // 获取xy位置  
	void set_ellipse();			// 求误差椭圆	
	vector<vector<double>> get_A();	//获取A矩阵
	vector<double> get_b(); //获取b矩阵
	double set_z();				// 根据定位点获得z的估计值 
	double func_fan(double sig11, double sig12, double sig22);	//根据x轴 y轴 方差得到偏角.


public:
	PLocation();			//无参构造函数
	PLocation(vector<DFdata> dflist_in);	//有参构造函数
	void input(vector<double> x_in, vector<double> y_in, vector<double> z_in, vector<double> alpha_in);		//从原始数据输入，与output成对
	vector<vector<double>> output();
    vector<vector<double>> outputByKmeans();																													//output 与input成对
	void set_dflist(vector<DFdata> dflist_in); // 更新dflist并初始化
	void set_Pe(double pe_in);	//设置置信度
	double get_Pe();			//获取当前置信度
	int get_n();			//获取 当前运算数据的数量
	vector<double> get_locxyz(); // 获取坐标xyz		  
	vector<double> get_ellipse(); //获取误差椭圆	
	vector<DFdata> get_dflist();	//获取dflist
	bool isValid();				//是否有效
	double AverageDistance;		//平均距离
};
#endif