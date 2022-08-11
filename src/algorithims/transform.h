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


// �Ƕ� ���� ת��
double hu2angle(double hu);
double angle2hu(double angle);
double to360(double angle);

// �󽻵�  ��λ�ǽǶ�
vector<double> cross_point(double x1, double y1, double alpha1,
                           double x2, double y2, double alpha2);

// �ж��Ƿ����������Զ
bool is_line_point(double platx, double platy, double alpha, double tarx, double tary);

//������֮��Ƕ�
double get_angle(double platx, double platy, double tarx, double tary);

//time
struct tm create_time(int year, int month, int date, int hour, int min, int sec);
string to_str_time(struct tm t);
//˵�����������t_early �� t_late�磬������ֵ�����򷵻ظ�ֵ
long sub_time(struct tm t_early, struct tm t_late);

bool same_point(double x1, double y1, double x2, double y2);


//alpha1 ����Ϊ����
double getPointLineDistance(double x1, double y1, double alpha1, double x2, double y2);


//Map �ӿ� Ƕ��ƽ̨ʱע�͵�
double getTwoPointDistance(double x1, double y1, double x2, double y2);


#endif //TRANSFORM_H
