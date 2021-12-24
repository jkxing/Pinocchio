#include "pinocchioApi.h"
#include "deriv.h"

struct RP //information for refined embedding
{
    RP(TreeType* inD, const Skeleton& inSk, const vector<Vector3>& medialSurface)
        : distanceField(inD), given(inSk)
    {
        vector<Vec3Object> mpts;
        for (int i = 0; i < (int)medialSurface.size(); ++i)
            mpts.push_back(medialSurface[i]);

        medProjector = ObjectProjector<3, Vec3Object>(mpts);
    }

    TreeType* distanceField;
    const Skeleton& given;
    ObjectProjector<3, Vec3Object> medProjector;
};

template<class Real> Real computeFineError(const vector<Vector<Real, 3> >& match, RP* rp) {
    Real out = Real();
    for (auto i = 1; i < match.size(); ++i) {
        int prev = rp->given.fPrev()[i];
        Real surfPenalty = Real();
        Real lenPenalty = Real();
        Real anglePenalty = Real();
        Real symPenalty = Real();
        const int samples = 10;
        for (int k = 0; k < samples; ++k) {
            double frac = double(k) / double(samples);
            Vector<Real, 3> cur = match[i] * Real(1. - frac) + match[prev] * Real(frac);
            Vector3 m = rp->medProjector.project(cur);
            Real medDist = (cur - Vector<Real, 3>(m)).length();
            Real surfDist = -rp->distanceField->locate(cur)->evaluate(cur);
            Real penalty = SQR(min(medDist, Real(0.001) + max(Real(0.), Real(0.05) - surfDist)));
            if (penalty > Real(SQR(0.003)))
                surfPenalty += Real(1. / double(samples)) * penalty;
        }

        //---------------length
        Real optDistSq = (rp->given.fGraph().verts[i] - rp->given.fGraph().verts[prev]).lengthsq();
        Real distSq = SQR(max(Real(-10.), (match[i] - match[prev]) *
            (rp->given.fGraph().verts[i] - rp->given.fGraph().verts[prev]))) / optDistSq;
        lenPenalty = max(Real(.5), (Real(0.0001) + optDistSq) / (Real(0.0001) + distSq));

        //---------------sym
        if (rp->given.fSym()[i] != -1) {
            int s = rp->given.fSym()[i];
            int sp = rp->given.fPrev()[s];

            Real sDistSq = (match[s] - match[sp]).lengthsq();
            symPenalty = max(Real(1.05), max(distSq / (Real(0.001) + sDistSq), sDistSq / (Real(0.001) + distSq)));
        }

        //--------------angle
        if (distSq > Real(1e-16)) {
            Vector<Real, 3> curDir = (match[i] - match[prev]).normalize();
            Vector<Real, 3> skelDir = (rp->given.fGraph().verts[i] - rp->given.fGraph().verts[prev]).normalize();
            if (curDir * skelDir < Real(1. - 1e-8))
                anglePenalty = Real(0.5) * acos(curDir * skelDir);
            anglePenalty = CUBE(Real(0.3) + anglePenalty);
            if (curDir * skelDir < Real(0.))
                anglePenalty *= 10.;
        }

        out += Real(15000.) * surfPenalty + Real(0.25) * lenPenalty + Real(2.0) * anglePenalty + symPenalty;
    }
    return out;
}

vector<Vector3> optimizeEmbedding1D(vector<Vector3> fineEmbedding, vector<Vector3> dir, RP* rp) {
    double lr = 0.001;
    for (auto i = 0; i < fineEmbedding.size(); ++i) {
        lr += dir[i].lengthsq();
    }
    lr = 0.0005 / sqrt(lr);
    double prevErr = -1e10;
    int count = 0;
    while (++count) {
        double curErr = computeFineError(fineEmbedding, rp);
        if (prevErr==-1e10||curErr < prevErr) {
            lr *= 2.;
            for (auto i = 0; i < fineEmbedding.size(); ++i) {
                fineEmbedding[i] += dir[i] * lr;
            }
            prevErr = curErr;
        }
        else {
            if (count > 2) {
                for (auto i = 0; i < fineEmbedding.size(); ++i) {
                    fineEmbedding[i] -= dir[i] * lr;
                }
            }
            break;
        }
    }
    return fineEmbedding;
}

vector<Vector3> m_refineEmbedding(TreeType* distanceField, const vector<Vector3>& medialSurface, const vector<Vector3>& initialEmbedding, const Skeleton& skeleton)
{
    RP rp(distanceField, skeleton, medialSurface);
    vector<Vector3> result(initialEmbedding);
    int iter_time = 10;
    int num = initialEmbedding.size();
    typedef Deriv<double, 6> DType;
    typedef Deriv<double, -1> DType1;
    for (auto k = 0; k < iter_time; k++) {
        Debugging::out() << "E = " << computeFineError(result, &rp) << endl;
        for (int j = 0; j < 2; ++j) {
            vector<Vector<DType1, 3> > dMatch(num);
            for (auto i = 0; i < num; ++i) {
                dMatch[i][0] = DType1(result[i][0], i * 3);
                dMatch[i][1] = DType1(result[i][1], i * 3 + 1);
                dMatch[i][2] = DType1(result[i][2], i * 3 + 2);
            }
            DType1 err = computeFineError(dMatch, &rp) + (DType1() * DType1(0., 3 * num));
            vector<Vector3> dir(num);
            for (auto i = 0; i < num; ++i) {
                dir[i][0] = -err.getDeriv(i * 3);
                dir[i][1] = -err.getDeriv(i * 3 + 1);
                dir[i][2] = -err.getDeriv(i * 3 + 2);
            }
            result = optimizeEmbedding1D(result, dir, &rp);
        }
    }
    return result;
}