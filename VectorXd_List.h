#include <stdio.h>
#include <stdlib.h>
#include<fstream>
#include<typeinfo>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <time.h>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

typedef struct VectorXd_List
{
	VectorXd V;
	struct VectorXd_List* next;
}VL;

//����һ���ڵ� 
VL* create_node(VectorXd* Vdata)
{
	//��ÿ���ڵ����ṹ��һ���Ŀռ��С 
	VL* p = (struct VectorXd_List*)malloc(sizeof(VL));
	if (NULL == p)
	{
		printf("malloc error!\n");
		return NULL;
	}
	//���ڽṹ����δ��ʼ����ʱ��һ���������ݣ�����Ҫ�� 
	memset(p, 0, sizeof(VL));
	//��ʼ����һ���ڵ� 
	p->V = *Vdata;
	//���ڵ�ĺ��ָ������ΪNULL 
	p->next = NULL;

	return p;
}

//�����β�� 
void tail_insert(VL* pH, VL* new1)
{
	//��ȡ��ǰ��λ�� 
	VL* p = pH;
	//�����ǰλ�õ���һ���ڵ㲻Ϊ�� 
	while (NULL != p->next)
	{
		//�ƶ�����һ���ڵ� 
		p = p->next;
	}
	//�����������ѭ���������Ѿ�����NULL�����λ��
	//��ʱֱ�Ӱ��²���Ľڵ㸳ֵ��NULL���λ�� 
	p->next = new1;
}


void Fprint_node(const char* filename,VL* pH, int m)
{
	//��ȡ��ǰ��λ�� 
	VL* p = pH;
	//��ȡ��һ���ڵ��λ�� 
	p = p->next;
	FILE* fp;
	if ((fp = fopen(filename,"w")) == NULL) {
		printf("cannot open this file\n");
		exit(0);
	}

	//�����ǰλ�õ���һ���ڵ㲻Ϊ�� 
	while (NULL != p->next)
	{
		//(1)���ڵ�����д���ļ� 
		/*
		cout << "Vector:\n" << *(p->V) << endl;
		for (int i = 0; i < m; i++) {
			cout << "Vector("<<i<<"):\n" << (*(p->V))(i) << endl;
		}
		*/
		for (int i = 0; i < m; i++) {
			fprintf(fp, "%.11f\n", (p->V)(i));
		}
		//(2)�ƶ�����һ���ڵ�,���������Ϊ�棬���ظ�(1)����(2) 
		p = p->next;
	}
	//�����ǰλ�õ���һ���ڵ�Ϊ�գ����ӡ����
	//˵��ֻ��һ���ڵ� 
	//cout << "Final Vector:\n" << *(p->V) << endl;
	for (int i = 0; i < m; i++) {
		fprintf(fp, "%.11f\n", (p->V)(i));
	}
	fclose(fp);
}

//�����ͷ�� 
void top_insert(VL* pH, VL* new1)
{
	VL* p = pH;
	new1->next = p->next;
	p->next = new1;
}

//�����м����
/*
void mid_insert(VL* pH, VectorXd* Vdata, VL* new1)
{
	VL* p = pH;
	p = p->next;
	while (p->V != Vdata)
	{
		p = p->next;
	}
	new1->next = p->next;
	p->next = new1;

}
*/


//����ı��� 
void Print_node(VL* pH)
{
	//��ȡ��ǰ��λ�� 
	VL* p = pH;
	//��ȡ��һ���ڵ��λ�� 
	p = p->next;
	//�����ǰλ�õ���һ���ڵ㲻Ϊ�� 
	while (NULL != p->next)
	{
		//(1)��ӡ�ڵ������ 
		cout << "Vector:\n" << p->V<<endl;
		//(2)�ƶ�����һ���ڵ�,���������Ϊ�棬���ظ�(1)����(2) 
		p = p->next;
	}
	//�����ǰλ�õ���һ���ڵ�Ϊ�գ����ӡ����
	//˵��ֻ��һ���ڵ� 
	cout << "Final Vector:\n" << p->V<<endl;
}

//ɾ�������еĽڵ� 
int detele_list_node(VL* pH, VectorXd* Vdata)
{
	//��ȡ��ǰͷ�ڵ��λ�� 
	VL* p = pH;
	VL* prev = NULL;
	while (NULL != p->next)
	{
		//���浱ǰ�ڵ��ǰһ���ڵ��ָ�� 
		prev = p;
		//Ȼ���õ�ǰ��ָ����������ƶ� 
		p = p->next;
		//�жϣ��ҵ���Ҫɾ��������  
		if (p->V == *Vdata)
		{
			//���������һ������ͨ�ڵ㣬����һ����β�ڵ�
			if (p->next != NULL)  //��ͨ�ڵ����� 
			{
				prev->next = p->next;
				free(p);
			}
			else //β�ڵ����� 
			{
				prev->next = NULL; //�����β�ڵ����һ���ڵ��ָ����ָ��� 
				free(p);
			}
			return 0;
		}
	}
	printf("û��Ҫɾ���Ľڵ�\n");
	return -1;
}

void trave_list(VL* pH)
{
	//�����һ���ڵ��λ�� 
	VL* p = pH->next;
	VL* pBack;
	int i = 0;
	if (p->next == NULL || p == NULL)//�������ֻ��1���ڵ������
		return;

	while (NULL != p->next) //�������� 
	{
		//�����һ���ڵ����һ���ڵ� 
		pBack = p->next;
		//�ҵ���һ����Ч�ڵ�,��ʵ����ͷָ�����һ���ڵ� 
		if (p == pH->next)
		{
			//��һ����Ч�ڵ�������һ���ڵ㣬����Ҫָ��NULL 
			p->next = NULL;
		}
		else
		{

			p->next = pH->next; //β������ 
		}
		pH->next = p; //Ϊ�˵ڶ����ڵ㰴��һ���ڵ�һ���������Դ����� 
		p = pBack; //����һ���ڵ� 
	}
	top_insert(pH, p); //�������һ���ڵ� 
}
