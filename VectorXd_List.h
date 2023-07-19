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

//创建一个节点 
VL* create_node(VectorXd* Vdata)
{
	//给每个节点分配结构体一样的空间大小 
	VL* p = (struct VectorXd_List*)malloc(sizeof(VL));
	if (NULL == p)
	{
		printf("malloc error!\n");
		return NULL;
	}
	//由于结构体在未初始化的时候一样是脏数据，所以要清 
	memset(p, 0, sizeof(VL));
	//初始化第一个节点 
	p->V = *Vdata;
	//将节点的后继指针设置为NULL 
	p->next = NULL;

	return p;
}

//链表的尾插 
void tail_insert(VL* pH, VL* new1)
{
	//获取当前的位置 
	VL* p = pH;
	//如果当前位置的下一个节点不为空 
	while (NULL != p->next)
	{
		//移动到下一个节点 
		p = p->next;
	}
	//如果跳出以上循环，所以已经到了NULL的这个位置
	//此时直接把新插入的节点赋值给NULL这个位置 
	p->next = new1;
}


void Fprint_node(const char* filename,VL* pH, int m)
{
	//获取当前的位置 
	VL* p = pH;
	//获取第一个节点的位置 
	p = p->next;
	FILE* fp;
	if ((fp = fopen(filename,"w")) == NULL) {
		printf("cannot open this file\n");
		exit(0);
	}

	//如果当前位置的下一个节点不为空 
	while (NULL != p->next)
	{
		//(1)将节点数据写入文件 
		/*
		cout << "Vector:\n" << *(p->V) << endl;
		for (int i = 0; i < m; i++) {
			cout << "Vector("<<i<<"):\n" << (*(p->V))(i) << endl;
		}
		*/
		for (int i = 0; i < m; i++) {
			fprintf(fp, "%.11f\n", (p->V)(i));
		}
		//(2)移动到下一个节点,如果条件仍为真，则重复(1)，再(2) 
		p = p->next;
	}
	//如果当前位置的下一个节点为空，则打印数据
	//说明只有一个节点 
	//cout << "Final Vector:\n" << *(p->V) << endl;
	for (int i = 0; i < m; i++) {
		fprintf(fp, "%.11f\n", (p->V)(i));
	}
	fclose(fp);
}

//链表的头插 
void top_insert(VL* pH, VL* new1)
{
	VL* p = pH;
	new1->next = p->next;
	p->next = new1;
}

//链表中间插入
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


//链表的遍历 
void Print_node(VL* pH)
{
	//获取当前的位置 
	VL* p = pH;
	//获取第一个节点的位置 
	p = p->next;
	//如果当前位置的下一个节点不为空 
	while (NULL != p->next)
	{
		//(1)打印节点的数据 
		cout << "Vector:\n" << p->V<<endl;
		//(2)移动到下一个节点,如果条件仍为真，则重复(1)，再(2) 
		p = p->next;
	}
	//如果当前位置的下一个节点为空，则打印数据
	//说明只有一个节点 
	cout << "Final Vector:\n" << p->V<<endl;
}

//删除链表中的节点 
int detele_list_node(VL* pH, VectorXd* Vdata)
{
	//获取当前头节点的位置 
	VL* p = pH;
	VL* prev = NULL;
	while (NULL != p->next)
	{
		//保存当前节点的前一个节点的指针 
		prev = p;
		//然后让当前的指针继续往后移动 
		p = p->next;
		//判断，找到了要删除的数据  
		if (p->V == *Vdata)
		{
			//两种情况，一种是普通节点，还有一种是尾节点
			if (p->next != NULL)  //普通节点的情况 
			{
				prev->next = p->next;
				free(p);
			}
			else //尾节点的情况 
			{
				prev->next = NULL; //将这个尾节点的上一个节点的指针域指向空 
				free(p);
			}
			return 0;
		}
	}
	printf("没有要删除的节点\n");
	return -1;
}

void trave_list(VL* pH)
{
	//保存第一个节点的位置 
	VL* p = pH->next;
	VL* pBack;
	int i = 0;
	if (p->next == NULL || p == NULL)//空链表或只有1个节点的链表
		return;

	while (NULL != p->next) //遍历链表 
	{
		//保存第一个节点的下一个节点 
		pBack = p->next;
		//找到第一个有效节点,其实就是头指针的下一个节点 
		if (p == pH->next)
		{
			//第一个有效节点就是最后一个节点，所以要指向NULL 
			p->next = NULL;
		}
		else
		{

			p->next = pH->next; //尾部连接 
		}
		pH->next = p; //为了第二个节点按第一个节点一样操作，以此类推 
		p = pBack; //走下一个节点 
	}
	top_insert(pH, p); //插入最后一个节点 
}
