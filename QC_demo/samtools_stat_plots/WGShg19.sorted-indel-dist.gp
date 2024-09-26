
        set terminal png size 600,400 truecolor
        set output "./output/QC/samtools_stat_plots/WGShg19.sorted-indel-dist.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style increment user
        set ylabel "Indel count [log]"
        set xlabel "Indel length"
        set y2label "Insertions/Deletions ratio"
        set log y
        set y2tics nomirror
        set ytics nomirror
        set title "WGShg19.sorted.sortedbamfilename.stats" noenhanced
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	49288
2	10343
3	4086
4	4019
5	4848
6	520
7	196
8	201
9	112
10	107
11	34
12	21
13	5
14	0
15	0
16	0
17	0
18	0
19	0
20	0
end
1	63833
2	14509
3	5160
4	5344
5	2018
6	949
7	354
8	439
9	197
10	345
11	137
12	113
13	63
14	55
15	31
16	23
17	6
18	8
19	2
20	2
end
1	0.772140
2	0.712868
3	0.791860
4	0.752058
5	2.402379
6	0.547945
7	0.553672
8	0.457859
9	0.568528
10	0.310145
11	0.248175
12	0.185841
13	0.079365
14	0.000000
15	0.000000
16	0.000000
17	0.000000
18	0.000000
19	0.000000
20	0.000000
end
