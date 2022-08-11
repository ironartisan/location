//
// Created by ChenYL on 2022/8/3.
//

#ifndef TRANSFORM_H
#define TRANSFORM_H
#include <cmath>
#include <vector>
#include <cstdio>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#define PI acos(-1)

using namespace std;


// 角度 弧度 转换
double hu2angle(double hu);
double angle2hu(double angle);
double to360(double angle);

// 求交点  单位是角度
vector<double> cross_point(double x1, double y1, double alpha1,
                           double x2, double y2, double alpha2);

// 判断是否点与射线相差不远
bool is_line_point(double platx, double platy, double alpha, double tarx, double tary);

//求两点之间角度
double get_angle(double platx, double platy, double tarx, double tary);

//time
struct tm create_time(int year, int month, int date, int hour, int min, int sec);
string to_str_time(struct tm t);
//说明：如果符合t_early 比 t_late早，返回正值，否则返回负值
long sub_time(struct tm t_early, struct tm t_late);

bool same_point(double x1, double y1, double x2, double y2);


//alpha1 度数为输入
double getPointLineDistance(double x1, double y1, double alpha1, double x2, double y2);


//Map 接口 嵌入平台时注释掉
double getTwoPointDistance(double x1, double y1, double x2, double y2);


#endif //TRANSFORM_H
