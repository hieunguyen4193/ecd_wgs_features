
            set terminal png size 600,400 truecolor
            set output "./output/QC/WGShg19.sorted/samtools_stat_plots/WGShg19.sorted-acgt-cycles.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set style line 1 linecolor rgb "green"
            set style line 2 linecolor rgb "red"
            set style line 3 linecolor rgb "black"
            set style line 4 linecolor rgb "blue"
            set style increment user
            set ylabel "Base content [%]"
            set xlabel "Read Cycle"
            set yrange [0:100]
            set title "WGShg19.sorted.sortedbamfilename.stats" noenhanced
            plot '-' w l ti 'A', '-' w l ti 'C', '-' w l ti 'G', '-' w l ti 'T'
        2	20.80
3	23.67
4	31.75
5	29.96
6	30.60
7	30.73
8	32.13
9	32.58
10	30.84
11	29.31
12	28.98
13	29.05
14	29.40
15	29.55
16	29.66
17	29.49
18	29.33
19	29.31
20	29.28
21	29.18
22	29.07
23	29.32
24	29.48
25	29.52
26	29.36
27	29.27
28	29.19
29	28.97
30	29.05
31	29.22
32	29.38
33	29.47
34	29.44
35	29.40
36	29.30
37	29.33
38	29.42
39	29.35
40	29.29
41	29.17
42	29.01
43	29.27
44	29.45
45	29.67
46	29.68
47	29.49
48	29.47
49	29.39
50	29.29
51	29.19
end
2	29.16
3	26.24
4	18.13
5	19.96
6	19.37
7	19.23
8	17.82
9	17.34
10	19.10
11	20.68
12	20.97
13	20.85
14	20.52
15	20.45
16	20.31
17	20.47
18	20.59
19	20.59
20	20.66
21	20.78
22	20.89
23	20.64
24	20.46
25	20.41
26	20.61
27	20.71
28	20.75
29	20.94
30	20.91
31	20.74
32	20.60
33	20.45
34	20.48
35	20.56
36	20.67
37	20.65
38	20.50
39	20.56
40	20.63
41	20.78
42	20.93
43	20.65
44	20.49
45	20.29
46	20.35
47	20.47
48	20.45
49	20.55
50	20.66
51	20.81
end
2	29.14
3	26.29
4	18.12
5	19.93
6	19.30
7	19.17
8	17.78
9	17.32
10	19.10
11	20.64
12	20.95
13	20.90
14	20.55
15	20.42
16	20.24
17	20.41
18	20.60
19	20.63
20	20.66
21	20.76
22	20.84
23	20.60
24	20.46
25	20.40
26	20.55
27	20.60
28	20.72
29	20.97
30	20.92
31	20.70
32	20.55
33	20.44
34	20.46
35	20.54
36	20.59
37	20.58
38	20.52
39	20.59
40	20.59
41	20.72
42	20.88
43	20.62
44	20.44
45	20.27
46	20.23
47	20.41
48	20.42
49	20.56
50	20.63
51	20.72
end
2	20.90
3	23.80
4	32.00
5	30.15
6	30.73
7	30.87
8	32.27
9	32.76
10	30.95
11	29.38
12	29.10
13	29.21
14	29.52
15	29.58
16	29.79
17	29.63
18	29.48
19	29.47
20	29.41
21	29.28
22	29.19
23	29.44
24	29.60
25	29.67
26	29.48
27	29.42
28	29.34
29	29.12
30	29.12
31	29.34
32	29.48
33	29.64
34	29.62
35	29.50
36	29.43
37	29.45
38	29.56
39	29.50
40	29.49
41	29.33
42	29.17
43	29.46
44	29.62
45	29.77
46	29.74
47	29.63
48	29.66
49	29.50
50	29.41
51	29.28
end