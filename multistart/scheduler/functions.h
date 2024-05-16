/*Pantelis Flouris, 1093507*/
/*Aggelos Menegatos, 1093426*/
/*Chrysafis Koltsakis, 1084671*/
/*Giorgos Amaxopoulos, 1093311*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <stdbool.h>
#include <time.h>
#include <sys/mman.h>
#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#define MAX 256


typedef enum status {
	NEW, RUNNING, STOPPED ,IO , EXITED
} status;

typedef struct Node {
	pid_t pid;
	char value[100];
	status status;
	struct Node* next;
} Node;

typedef struct Queue {
	Node *head, *tail;
} Queue;

Node* createNode(char *data, pid_t pid, int status);
Queue* createQueue(void);
bool Enqueue(Queue *queue, char* data, pid_t pid, int status);
Node* Dequeue(Queue* queue);

bool isNull(Queue* q);
double get_wtime(void);
#endif  // FUNCTIONS_H



double get_wtime(void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + (double)t.tv_usec * 1.0e-6;
}

struct Node* createNode(char *data, pid_t pid, int status){
    Node *node = mmap(NULL, sizeof(Node), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
	strncpy(node->value, (char*)data, sizeof(node->value) - 1);
	node->value[sizeof(node->value) - 1] = '\0';
	node->next=NULL;
	node->pid = pid;
	node->status = status;
	return node;
}

struct Queue* createQueue(){
    Queue *q = mmap(NULL, sizeof(Queue), PROT_READ 
					| PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
	q->head=NULL;
	q->tail=NULL;
	return q;
}

bool isNull(Queue* q){
	if(q->tail==NULL && q->head == NULL)
		return true;
	else return false;
}


bool Enqueue(Queue *queue, char* data,pid_t pid,int status){
	//create node
	Node* new = createNode(data, pid, status);
	if(new==NULL) return false;

	//check if queue is null, and add new node
	if(isNull(queue)){
		queue->tail = new;
		queue->head = new;
	}
	//else connect the new node as the last in the queue
	else {
	queue->tail->next= new;
	queue->tail = new;
	}
	return true;
}

Node* Dequeue(Queue* queue) {
    // Check if the queue is empty
    if (isNull(queue)) {
        return NULL;
    }

    // Save the pointer to the old head node
    Node* temp = queue->head;

    // Point to the new head node
    queue->head = queue->head->next;

    // If the queue is empty after dequeue, don't leave tail hanging
    if (queue->head == NULL) {
        queue->tail = NULL;
    }
    // Return the old head node
    return temp;
}

