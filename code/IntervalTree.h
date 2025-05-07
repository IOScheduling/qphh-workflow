#ifndef INTERVAL_TREE_H
#define INTERVAL_TREE_H

#include <vector>
#include <stack>
#include <algorithm>
using namespace std;

// �����������ڵ�
struct Interval {
    double startTime;
    double endTime;
    int taskId;

    Interval(double start, double end, int id) : startTime(start), endTime(end), taskId(id) {}
};

// ��������
class IntervalTree {
private:
    struct TreeNode {
        double startTime;
        double maxEndTime;
        vector<Interval> intervals; // �洢�����ڴ˽ڵ㷶Χ�ڵ�����
        TreeNode* left;
        TreeNode* right;

        TreeNode(double start) : startTime(start), maxEndTime(start), left(nullptr), right(nullptr) {}
    };

    TreeNode* root;

    // �ڲ����뺯����ʹ�õ�����ʽ
    void insertImpl(double startTime, double endTime, int taskId) {
        TreeNode* node = root;
        TreeNode* parent = nullptr;

        while (node != nullptr) {
            parent = node;
            node->maxEndTime = max(node->maxEndTime, endTime);
            if (startTime < node->startTime) {
                node = node->left;
            }
            else {
                node = node->right;
            }
        }

        // �����µĽڵ�
        TreeNode* newNode = new TreeNode(startTime);
        newNode->maxEndTime = endTime;
        newNode->intervals.push_back(Interval(startTime, endTime, taskId));

        if (parent == nullptr) {
            root = newNode;
        }
        else if (startTime < parent->startTime) {
            parent->left = newNode;
        }
        else {
            parent->right = newNode;
        }
    }

    // ��ѯ����������ص������䣬ʹ�õ�����ʽ
    vector<int> queryImpl(double startTime, double endTime) {
        vector<int> result;
        stack<TreeNode*> nodeStack;
        TreeNode* node = root;
        while (node != nullptr || !nodeStack.empty()) {
            if (node != nullptr) {
                if (node->maxEndTime >= startTime) {
                    // �����ǰ����Ľ���ʱ����ڲ�ѯ����Ŀ�ʼʱ�䣬���п������ص�
                    for (const auto& interval : node->intervals) {
                        if (!(interval.endTime < startTime || interval.startTime > endTime)) {
                            result.push_back(interval.taskId);
                        }
                    }
                    nodeStack.push(node);
                    node = node->left;
                }
                else {
                    node = node->right;
                }
            }
            else {
                node = nodeStack.top();
                nodeStack.pop();
            }
        }
        return result;
    }

    // ɾ�����䣬ʹ�õ�����ʽ
    void removeImpl(double startTime, double endTime, int taskId) {
        stack<TreeNode*> nodeStack;
        TreeNode* node = root;
        TreeNode* parent = nullptr;

        while (node != nullptr) {
            if (node->startTime == startTime && node->maxEndTime == endTime) {
                auto& intervals = node->intervals;
                intervals.erase(remove_if(intervals.begin(), intervals.end(),
                    [startTime, endTime, taskId](const Interval& interval) {
                        return interval.startTime == startTime && interval.endTime == endTime && interval.taskId == taskId;
                    }), intervals.end());

                // ����������ʱ��
                node->maxEndTime = node->startTime;
                for (auto& interval : node->intervals) {
                    node->maxEndTime = max(node->maxEndTime, interval.endTime);
                }

                break;
            }
            parent = node;
            if (startTime < node->startTime) {
                node = node->left;
            }
            else {
                node = node->right;
            }
        }
    }

public:
    IntervalTree() : root(nullptr) {}

    // �����ӿڵ���˽��ʵ�ֺ���
    void insert(double startTime, double endTime, int taskId) {
        insertImpl(startTime, endTime, taskId);
    }

    vector<int> query(double startTime, double endTime) {
        return queryImpl(startTime, endTime);
    }

    void remove(double startTime, double endTime, int taskId) {
        removeImpl(startTime, endTime, taskId);
    }
};

#endif // INTERVAL_TREE_H
