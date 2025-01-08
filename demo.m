%% read testfile and groundtruth
data = readmatrix("./testdata/testfile2/test_file.csv");
gt1 = readmatrix("./testdata/testfile2/gt1.csv");
gt2 = readmatrix("./testdata/testfile2/gt2.csv");


%% run SCAN
[Omega, pred] = scan_wrapper(data);

%% Compute Jaccard similarities
accuracy_t1 = jaccard(pred{1}, gt1);
accuracy_t2 = jaccard(pred{2}, gt2);

%% visualize result
figure; imagesc(data)
figure; imagesc(pred{1});
figure; imagesc(pred{2});