#pragma GCC optimize(1)
#pragma GCC optimize(2)
#pragma GCC optimize(3)
#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#pragma GCC optimize("-fgcse")
#pragma GCC optimize("-fgcse-lm")
#pragma GCC optimize("-fipa-sra")
#pragma GCC optimize("-ftree-pre")
#pragma GCC optimize("-ftree-vrp")
#pragma GCC optimize("-fpeephole2")
#pragma GCC optimize("-ffast-math")
#pragma GCC optimize("-fsched-spec")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("-falign-jumps")
#pragma GCC optimize("-falign-loops")
#pragma GCC optimize("-falign-labels")
#pragma GCC optimize("-fdevirtualize")
#pragma GCC optimize("-fcaller-saves")
#pragma GCC optimize("-fcrossjumping")
#pragma GCC optimize("-fthread-jumps")
#pragma GCC optimize("-funroll-loops")
#pragma GCC optimize("-freorder-blocks")
#pragma GCC optimize("-fschedule-insns")
#pragma GCC optimize("inline-functions")
#pragma GCC optimize("-ftree-tail-merge")
#pragma GCC optimize("-fschedule-insns2")
#pragma GCC optimize("-fstrict-aliasing")
#pragma GCC optimize("-falign-functions")
#pragma GCC optimize("-fcse-follow-jumps")
#pragma GCC optimize("-fsched-interblock")
#pragma GCC optimize("-fpartial-inlining")
#pragma GCC optimize("no-stack-protector")
#pragma GCC optimize("-freorder-functions")
#pragma GCC optimize("-findirect-inlining")
#pragma GCC optimize("-fhoist-adjacent-loads")
#pragma GCC optimize("-frerun-cse-after-loop")
#pragma GCC optimize("inline-small-functions")
#pragma GCC optimize("-finline-small-functions")
#pragma GCC optimize("-ftree-switch-conversion")
#pragma GCC optimize("-foptimize-sibling-calls")
#pragma GCC optimize("-fexpensive-optimizations")
#pragma GCC optimize("inline-functions-called-once")
#pragma GCC optimize("-fdelete-null-pointer-checks")
#include<iostream>
#include<fstream>
#include<cstring>
#include<string>
#include<vector>
#include<cstdio>
#include<unordered_map>
#include<unordered_set>
#include<set>
#include<map>
#include<algorithm>
#include<chrono>
#include<array>
#include<queue>
#include<ctime>
#include<stdlib.h>
#include<sstream>
using namespace std;

template<typename T>
void hash_combine(size_t& seed, T const& v)
{
	seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
struct pair_hash
{
	template<typename T1, typename T2>
	size_t operator()(std::pair<T1, T2>const& p) const
	{
		size_t seed = 0;
		hash_combine(seed, p.first);
		hash_combine(seed, p.second);
		return seed;
	}
};

typedef unsigned int UI;
typedef long long LL;
typedef pair<int, int> PII;
typedef pair<UI, int>PUI;
typedef pair<int, pair<int, int> > PIII;
#define yx second.first
#define yy second.second
#define x first
#define y second

const UI INF = 0xffffffffll;
int vern = 0, kmax = 0;
int m = 0;
UI tmax = 0;
int tmin = INT_MAX;

vector<int> core, recCore, d;
vector<vector<PUI>> H_edges;
struct Edge {
	int u, v;
	LL t;
	bool operator<(const Edge& E)const
	{
		if (t != E.t)return t < E.t;
		if (u != E.u)return u < E.u;
		return v < E.v;
	}
};
typedef pair<UI, UI>PUU;
vector<Edge>edges;
struct Query {
	int st, ed, k;
};
vector<Query>q;
int qcnt = 0;
vector<int> vers;
vector<vector<vector<PUU>>> PHC_Index;
vector<vector<UI>>CT;
vector<unordered_map<int, int>>CN;
vector<vector<unordered_map<int, int>>>CTN;
unordered_map<UI, unordered_set<PII, pair_hash>>recedges;
vector<bool> inqueue;
bool CTN_Used[5000010];
vector<unordered_map<PUU, int, pair_hash>>recMinK;
vector<unordered_map<PUU, int, pair_hash>>recMaxK;
vector<unordered_map<PUU, int, pair_hash>>recId;
struct Node
{
	int u, te, kmin, kmax;
	vector<vector<PII>>ne;
};
vector<Node>nodes;
void delmp(unordered_map<int, int>& rec)
{
	unordered_map<int, int> tmp;
	rec.swap(tmp);
}
int get_id(int u)
{
	return (lower_bound(vers.begin(), vers.end(), u) - vers.begin()) + 1;
}

void loadgraph(const char* name)
{
	ifstream fin(name, ios::in);
	if (fin.is_open() == false) { printf("open graph %s fail\n", name); exit(1); }
	unordered_set<int>rectime;
	int msize = 0;
	string l;
	stringstream sstream;
	LL ttmax = 0;
	while (getline(fin, l))
	{
		sstream.clear();
		sstream.str(l);
		msize++;
		int u, v, t;
		sstream >> u >> v >> t;
		vers.push_back(u), vers.push_back(v);
		if (u == v)continue;
		rectime.insert(t);
		tmin = min(tmin, t);
		if (u > v)swap(u, v);
		edges.push_back({ u, v, t });
	}
	m = edges.size();
	sort(edges.begin(), edges.end());
	int k = 0;
	for (int i = 1; i < m; i++)
	{
		if (edges[i].u != edges[i - 1].u || edges[i].v != edges[i - 1].v || edges[i].t != edges[i - 1].t)
		{
			k++;
			edges[k] = edges[i];
		}
	}
	edges.resize(k + 1);
	m = k + 1;
	sort(vers.begin(), vers.end());
	vers.erase(unique(vers.begin(), vers.end()), vers.end());
	vern = vers.size();
	H_edges.resize(vern + 1);
	d.resize(vern + 1, 0);
	inqueue.resize(vern + 1, false);
	for (int i = 0; i < m; ++i)
	{
		ttmax = max(ttmax, edges[i].t + 1ll - tmin);
		edges[i].t -= tmin - 1;
		edges[i].u = get_id(edges[i].u), edges[i].v = get_id(edges[i].v);
		if (edges[i].u > edges[i].v) swap(edges[i].u, edges[i].v);
		tmax = max(tmax, (UI)edges[i].t);
		H_edges[edges[i].u].push_back({ edges[i].t, edges[i].v });
		H_edges[edges[i].v].push_back({ edges[i].t, edges[i].u });
	}
	for (int i = 0; i < m; i++)
	{
		recedges[edges[i].t].insert({ edges[i].u, edges[i].v });
	}
	printf("vern size: %d\nedge size: %d \ndifferent t is %d\n", vern, msize, (int)rectime.size());
	cout << "tmin : " << tmin << "\ntmax : " << tmax << "\nttmax : " << ttmax << "\n";
}


void coredecomp()
{
	core.resize(vern + 1, 0);
	int dmax = 0;
	unordered_set<PII, pair_hash> recEdge;
	for (int i = 0; i < m; i++)
	{
		int u = edges[i].u, v = edges[i].v;
		if (recEdge.count({ u,v }))
		{
			continue;
		}
		recEdge.insert({ u,v });
		d[u]++;
		d[v]++;
		dmax = max(dmax, d[u]);
		dmax = max(dmax, d[v]);
	}
	unordered_map<int, unordered_set<int>> recv;
	for (int u = 1; u <= vern; u++)
	{
		recv[d[u]].insert(u);
	}

	vector<int>state(vern + 1, 0);
	for (int k = 1; k <= dmax; k++)
	{
		if (!recv.count(k))continue;
		queue<int>q;
		for (int u : recv[k])
		{
			q.push(u);
		}
		while (!q.empty())
		{

			int u = q.front();
			q.pop();
			core[u] = k;
			for (PUI edge : H_edges[u])
			{
				int x = edge.y;
				if (state[x] == u)
					continue;
				state[x] = u;
				if (d[x] > k)
				{
					recv[d[x]].erase(x);
					d[x]--;
					if (d[x] > k)
					{
						recv[d[x]].insert(x);
					}
					else
					{
						q.push(x);
					}

				}
			}
		}
	}
	recCore = core;
	PHC_Index.resize(vern + 1);
	for (int i = 1; i <= vern; i++)
	{
		PHC_Index[i].resize(core[i]);
	}
	for (int i = 1; i <= vern; i++)
	{
		kmax = max(kmax, core[i]);
	}
	printf("kmax = %d\n", kmax);
}
void initCN()
{
	CN.resize(vern + 1);
	for (int i = 0; i < m; i++)
	{
		int u = edges[i].u, v = edges[i].v;
		UI t = edges[i].t;
		if (core[u] <= core[v])
		{
			if (!CN[u].count(v))
			{
				CN[u][v] = 0;
			}
			CN[u][v]++;
		}
		if (core[v] <= core[u])
		{
			if (!CN[v].count(u))
			{
				CN[v][u] = 0;
			}
			CN[v][u]++;
		}
	}
}
void clearCN()
{
	for(int i = 1;i <= vern;i++)
	{
		delmp(CN[i]);
	}
	vector<unordered_map<int, int>>tCN;
	CN.swap(tCN);
}
int localCore(int u, UI ts, UI te)
{
	unordered_set<int>rec;
	vector<int>cnt(core[u] + 1, 0);
	for (PUI edge : H_edges[u])
	{
		int x = edge.y;
		if (edge.x > te)break;
		if (edge.x >= ts)
		{
			if (rec.count(x))continue;
			rec.insert(x);
			cnt[min(core[u], core[x])]++;
		}
	}
	int cd = 0;
	for (int k = core[u]; k >= 1; k--)
	{
		cd = cd + cnt[k];
		if (cd >= k)return k;
	}
	return 0;
}
void CalCN(int u, UI ts, UI te)
{
	delmp(CN[u]);
	for (PUI edge : H_edges[u])
	{
		int x = edge.y;
		if (edge.x > te)break;
		if (edge.x >= ts && core[x] >= core[u])
		{
			CN[u][x]++;
		}
	}
}
void deledges(UI t, UI ts, UI te)
{
	if (!recedges.count(t))return;
	queue<int>q;
	for (PUI edge : recedges[t])
	{
		int u = edge.x, v = edge.y;
		if (core[u] <= core[v])
		{
			CN[u][v]--;
			if (CN[u][v] == 0)
			{
				CN[u].erase(v);
				if (CN[u].size() < core[u])
				{
					if (!inqueue[u])
					{
						inqueue[u] = true;
						q.push(u);

					}
				}
			}
		}
		if (core[v] <= core[u])
		{
			CN[v][u]--;
			if (CN[v][u] == 0)
			{
				CN[v].erase(u);
				if (CN[v].size() < core[v])
				{
					if (!inqueue[v])
					{
						inqueue[v] = true;
						q.push(v);
					}
				}
			}
		}
	}
	while (!q.empty())
	{
		int u = q.front();
		q.pop();
		inqueue[u] = false;
		int oldCore = core[u];
		core[u] = localCore(u, ts, te);
		CalCN(u, ts, te);
		for (int i = oldCore; i > core[u]; i--)
		{
			PHC_Index[u][i - 1].push_back({ ts,te + 1 });
		}
		for (PUI edge : H_edges[u])
		{
			if (edge.x > te)
				break;
			int x = edge.y;
			if (edge.x >= ts)
			{
				if (core[x] <= oldCore && core[x] > core[u])
				{
					if (CN[x].count(u))
					{
						CN[x].erase(u);
						if (CN[x].size() < core[x])
						{
							if (!inqueue[x])
							{
								inqueue[x] = true;
								q.push(x);
							}
						}
					}
				}
			}
		}
	}
}
void initCTN()
{
	core = recCore;
	CTN.resize(vern + 1);
	CT.resize(vern + 1);
	for (int u = 1; u <= vern; u++)
	{
		CT[u].resize(recCore[u]);
		for (int k = 1; k <= core[u]; k++)
		{
			if (PHC_Index[u][k - 1].size() == 0)
			{
				printf("error in %d %d\n", u, k);
			}
			CT[u][k - 1] = PHC_Index[u][k - 1][0].y;
		}
	}
	for (int u = 1; u <= vern; u++)
	{
		CTN[u].resize(core[u]);
		for (int k = 1; k <= core[u]; k++)
		{
			UI bound = PHC_Index[u][k - 1][0].y;
			for (PUI edge : H_edges[u])
			{
				int x = edge.y;
				if (edge.x > bound)break;
				if (core[x] >= k && PHC_Index[x][k - 1][0].y <= bound)
				{
					CTN[u][k - 1][x]++;
				}
			}
		}
	}
}
void initCTN(int u)
{
	return;
	
}
void DelCTN(int u, int v, int k, UI t, unordered_set<PII, pair_hash>& q)
{
	if (max(CT[v][k - 1], t) <= CT[u][k - 1])
	{
		CTN[u][k - 1][v]--;
		if (CTN[u][k - 1][v] == 0)
		{
			CTN[u][k - 1].erase(v);
			if (CTN[u][k - 1].size() < k)
			{
				q.insert({ u,k });
			}
		}
	}
}
UI LocalCT(int u, UI ts, int k)
{
	priority_queue<UI>q;
	unordered_set<int>rec;
	vector<PUI>& u_edge = H_edges[u];
	int l = 0, r = u_edge.size();
	while (l < r)
	{
		int mid = l + r >> 1;
		if (u_edge[mid].x >= ts)
			r = mid;
		else l = mid + 1;
	}
	for (int i = l; i < u_edge.size(); i++)
	{
		UI t = u_edge[i].x;
		int v = u_edge[i].y;
		if (q.size() == k && t >= q.top())break;
		if (rec.count(v))continue;
		rec.insert(v);
		if (core[v] < k)continue;
		UI t_v = max(t, CT[v][k - 1]);
		if (q.size() == k)
		{
			if (q.top() > t_v)
			{
				q.pop();
				q.push(t_v);
			}
		}
		else
		{
			q.push(t_v);
		}
	}
	if (q.size() == k)return q.top();
	int t = min(core[u], k - 1);
	core[u] = min(core[u], k - 1);
	return INF;

}
void calCTN(int u, int k, UI ts)
{
	if (core[u] < k)return;
	delmp(CTN[u][k - 1]);
	vector<PUI>& u_edge = H_edges[u];
	int l = 0, r = u_edge.size();
	while (l < r)
	{
		int mid = l + r >> 1;
		if (u_edge[mid].x >= ts)
			r = mid;
		else l = mid + 1;
	}
	UI bound = CT[u][k - 1];
	for (int i = l; i < u_edge.size(); i++)
	{
		UI t = u_edge[i].x;
		int v = u_edge[i].y;
		if (t > bound)break;
		if (core[v] < k)continue;
		UI t_v = max(t, CT[v][k - 1]);
		if (t_v <= bound)
			CTN[u][k - 1][v]++;
	}
}
void buildPHC()
{
	auto tstart = clock();
	coredecomp();
	printf("finish coredecomp when %lf s\n", ((clock() - tstart) / (double)CLOCKS_PER_SEC));
	initCN();
	printf("finish initCN when %lf\n", ((clock() - tstart) / (double)CLOCKS_PER_SEC));
	for (UI t = tmax; t > 0; t--)
	{
		deledges(t, 1, t - 1);
	}
	printf("finish deledges when %lf s\n", ((clock() - tstart) / (double)CLOCKS_PER_SEC));
	clearCN();
	initCTN();
	printf("finish initCTN when %lf s\n", ((clock() - tstart) / (double)CLOCKS_PER_SEC));
	unordered_set<PII, pair_hash> q;
	for (UI ts = 1; ts <= tmax; ts++)
	{
		if (!recedges.count(ts - 1))continue;
		for (PII edge : recedges[ts - 1])
		{
			int u = edge.x;
			int v = edge.y;
			initCTN(u);
			initCTN(v);
			for (int k = min(core[u], core[v]); k >= 1; k--)
			{
				DelCTN(u, v, k, ts - 1, q);
				DelCTN(v, u, k, ts - 1, q);
			}
		}
		while (!q.empty())
		{
			PII t = *(q.begin());
			q.erase(t);
			int u = t.x, k = t.y;
			UI oldCT = CT[u][k - 1];
			CT[u][k - 1] = LocalCT(u, ts, k);
			if (PHC_Index[u][k - 1][PHC_Index[u][k - 1].size() - 1].x != ts)
				PHC_Index[u][k - 1].push_back({ ts,CT[u][k - 1] });
			else PHC_Index[u][k - 1][PHC_Index[u][k - 1].size() - 1].y = CT[u][k - 1];
			calCTN(u, k, ts);
			unordered_set<int> rec_vertex;
			vector<PUI>& u_edge = H_edges[u];
			int l = 0, r = u_edge.size();
			while (l < r)
			{
				int mid = l + r >> 1;
				if (u_edge[mid].x >= ts)
					r = mid;
				else l = mid + 1;
			}
			UI bound = CT[u][k - 1];
			for (int i = l; i < u_edge.size(); i++)
			{
				UI t = u_edge[i].x;
				int v = u_edge[i].y;
				if (t >= bound)break;
				if (rec_vertex.count(v))continue;
				rec_vertex.insert(v);
				if (recCore[v] < k)continue;
				if (max(oldCT, t) <= CT[v][k - 1] && CT[v][k - 1] < CT[u][k - 1])
				{
					initCTN(v);
					CTN[v][k - 1].erase(u);
					if (CTN[v][k - 1].size() < k)
					{
						q.insert({ v,k });
					}
				}
			}

		}
	}
	printf("finish build PHC Index when %lf s\n", ((clock() - tstart) / (double)CLOCKS_PER_SEC));
}
void storePHC()
{
	//ios::in
	ofstream fout("phc_file.txt", ios::out);
	if (!fout.is_open())
	{
		printf("fail to open result file\n");
		return;
	}
	fout << vern << " " << tmax << " " << kmax << "\n";
	for (int i = 0; i < vern; i++)
	{
		if (i)
		{
			fout << " ";
		}
		fout << vers[i];
	}
	fout << "\n";
	for (int i = 1; i <= vern; i++)
	{
		fout << recCore[i] << "\n";
		for (int j = 0; j < recCore[i]; j++)
		{
			vector<PUU>& tmp = PHC_Index[i][j];
			int t = tmp.size();
			fout << t << "\n";
			for (int k = 0; k < t; k++)
			{
				fout << tmp[k].x << " " << tmp[k].y << "\n";
			}
		}
	}
	fout.close();
}
int main(int argc, char* argv[])
{
//	loadgraph("CollegeMsg.txt");
	loadgraph("SuperUser.txt");

	buildPHC();
	long long totalSize = 0;
	for (int i = 1; i <= vern; i++)
	{
		for (int j = 1; j < recCore[i]; j++)
			totalSize += PHC_Index[i][j].size() * 8;
	}
	int MB = 1 << 20; 
	printf("PHC_Index Size is : %.2lf MB\n", totalSize / (double)MB);
	storePHC();
	return 0;
}

