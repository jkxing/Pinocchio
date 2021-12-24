#include "pinocchioApi.h"
#include "deriv.h"
#include "debugging.h"

vector<Sphere> m_sampleMedialSurface(TreeType* distanceField, double tol)
{
	vector<Sphere> result;
	queue<OctTreeNode*> que;
	que.push(distanceField);
	auto getMinDot = [](TreeType* distanceField, Vector3 p, double r)
	{
		typedef Deriv<double, 3> D;
		typedef Vector<D, 3> VD;
		vector<Vector3> dfs;
		double mindot = 1.;
		for (int i = 0; i < 8; i++) {
			Vector3 f((i & 1) ? r : -r, (i >> 1 & 1) ? r : -r, (i >> 2 & 1) ? r : -r);
			f += p;
			VD vd = VD(D(f[0], 0), D(f[1], 1), D(f[2], 2));
			D res = distanceField->locate(f)->evaluate(vd);
			Vector3 df = Vector3(res.getDeriv(0), res.getDeriv(1), res.getDeriv(2)).normalize();
			for (auto& i : dfs) {
				mindot = min(mindot, df * i);
			}
			dfs.push_back(df);
		}
		return mindot;
	};
	while (!que.empty())
	{
		OctTreeNode* cur = que.front();
		que.pop();
		if (cur->getChild(0))
		{
			for (auto i = 0; i < 8; ++i)
				que.push(cur->getChild(i));
			continue;
		}

		Rect3 r = cur->getRect();
		double rad = r.getSize().length() / 2;
		Vector3 c = r.getCenter();
		double dot = getMinDot(distanceField, c, rad);
		if (dot > 0.)
			continue;

		double step = tol;
		vector<Vector3> points;
		double sz = r.getSize()[0];
		for (double x = 0; x <= sz; x += step) {
			for (double y = 0; y <= sz; y += step) {
				points.push_back(r.getLo() + Vector3(x, y, 0));
				if(y!=0) points.push_back(r.getLo() + Vector3(x, 0, y));
				if(x!=0&&y!=0) points.push_back(r.getLo() + Vector3(0, x, y));
			}
		}

		for (auto &p : points) {
			double dis = -distanceField->locate(p)->evaluate(p);
			if (dis <= 2. * step)
				continue;
			double dot = getMinDot(distanceField, p, step * 0.001);
			if (dot > 0.0)
				continue;
			result.push_back(Sphere(p, dis));
		}
	}
	sort(result.begin(), result.end(), [](const Sphere& x, const Sphere& y) {return x.radius > y.radius; });
	return result;
}

vector<Sphere> m_packSpheres(const vector<Sphere>& spheres, int maxSpheres) {
	
	vector<Sphere> result;
	for (auto &s:spheres) {
		bool flag = false;
		for (auto &r:result) {
			if ((r.center - s.center).length() < r.radius)
			{
				flag = true;
				break;
			}
		}
		if (flag) continue;
		result.push_back(s);
		if (result.size() >= maxSpheres)
			break;
	}
	return result;
}

PtGraph m_connectSamples(TreeType* distanceField, const vector<Sphere>& spheres) {
	PtGraph result;
	for (auto& s : spheres) {
		result.verts.push_back(s.center);
	}
	result.edges.resize(spheres.size());
	auto link = [&result](int i, int j) {
		result.edges[i].push_back(j);
		result.edges[j].push_back(i);
	};
	auto getMaxDist = [](TreeType* distanceField, Vector3 c1, Vector3 c2, double tol) {
		double mx = 1e-10;
		auto diff = (c2 - c1) / 100;
		for (auto i = 0; i < 101; ++i) {
			auto p = c1 + diff * i;
			mx = max(mx, distanceField->locate(p)->evaluate(p));
			if (mx > tol)
				break;
		}
		return mx;
	};
	for (auto i = 1; i < spheres.size();++i) {
		for (auto j = 0; j < i; ++j) {
			double dis = (spheres[i].center - spheres[j].center).length();
			if (dis < spheres[i].radius + spheres[j].radius)
			{
				link(i, j);
				continue;
			}
			double radsq = (spheres[i].center - spheres[j].center).lengthsq()/4;
			auto midP = (spheres[i].center + spheres[j].center) / 2;
			bool flag = false;
			for (auto k = 0; k < spheres.size(); ++k) {
				if (k == i || k == j) continue;
				if ((midP - spheres[k].center).lengthsq() < radsq)
				{
					flag = true;
					break;
				}
			}
			if (flag) continue;
			double maxAllowed = -0.5 * min(spheres[i].radius, spheres[j].radius);
			if (getMaxDist(distanceField, spheres[i].center, spheres[j].center, maxAllowed) < maxAllowed)
				link(i, j);
		}
	}
	return result;
}
