#include "pinocchioApi.h"
#include "debugging.h"

PinocchioOutput m_autorig(const Skeleton& given, const Mesh& m)
{
    int i;
    PinocchioOutput out;

    Mesh newMesh = prepareMesh(m);

    if (newMesh.vertices.size() == 0)
        return out;

    TreeType* distanceField = constructDistanceField(newMesh);

    vector<Sphere> medialSurface = sampleMedialSurface(distanceField);

    vector<Sphere> spheres = packSpheres(medialSurface);

    PtGraph graph = connectSamples(distanceField, spheres);

    vector<vector<int> > possibilities = computePossibilities(graph, spheres, given);

    vector<int> embeddingIndices = discreteEmbed(graph, spheres, given, possibilities);

    if (embeddingIndices.size() == 0) { //failure
        delete distanceField;
        return out;
    }

    vector<Vector3> discreteEmbedding = splitPaths(embeddingIndices, graph, given);

    vector<Vector3> medialCenters(medialSurface.size());
    for (i = 0; i < (int)medialSurface.size(); ++i)
        medialCenters[i] = medialSurface[i].center;

    out.embedding = refineEmbedding(distanceField, medialCenters, discreteEmbedding, given);

    VisTester<TreeType>* tester = new VisTester<TreeType>(distanceField);
    out.attachment = new Attachment(newMesh, given, out.embedding, tester);

    delete tester;
    delete distanceField;

    return out;
}

