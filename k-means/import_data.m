function [dataset, index, centroid, sizes] = ...
    import_data(numObjects,numAttributes,numClusters)
% function [dataset, index, centroid, sizes] = ...
%     import_data(numObjects,numAttributes,numClusters)
% Import data from binary files
% 'centroids.bin' 'dataset.bin' 'Index.bin' 'ClusterSize.bin'
% where
%   numObjects: elements of the initial dataset
%   numAttributes: number of attributes of the elements
%   numClusters: number of Clusters
%

fid = fopen('centroids.bin');
centroid = fread(fid, [numAttributes numClusters], 'double')';
fclose(fid);

fid = fopen('dataset.bin');
dataset = fread(fid, [numAttributes, numObjects], 'double')';
fclose(fid);

fid = fopen('Index.bin');
index = fread(fid, [numObjects 1], 'uint32');
fclose(fid);

fid = fopen('ClusterSize.bin');
sizes = fread(fid, [numClusters 1], 'uint32');
fclose(fid);


