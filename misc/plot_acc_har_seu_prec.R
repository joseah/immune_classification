library(tidyverse)
library(dsLib)


harmony <- data.table::fread(text = "	B cell	CD4+ T cell	CD8+ T cell	gd T cell	NKT cell	NK cell	C Monocyte	NC-Int Monocyte	cDC	pDC
B cell	327	0	0	0	0	0	3	0	0	0
CD4+ T cell	0	1505	9	3	0	0	0	0	1	0
CD8+ T cell	0	19	315	3	0	2	0	0	0	0
gd T cell	0	3	9	44	0	0	0	0	0	0
NKT cell	0	0	1	0	5	6	0	0	0	0
NK cell	0	0	1	0	1	315	0	0	0	0
C Monocyte	0	0	0	0	0	0	1106	32	5	0
NC-Int Monocyte	0	0	0	0	0	0	0	58	1	0
cDC	0	0	0	0	0	0	5	0	83	0
pDC	0	0	0	0	0	0	0	0	7	0", data.table = FALSE)


seurat <- data.table::fread(
"	B cell	CD4+ T cell	CD8+ T cell	gd T cell	NKT cell	NK cell	C Monocyte	NC-Int Monocyte	cDC	pDC
B cell	326	0	0	0	0	0	3	0	1	0
CD4+ T cell	0	1449	5	2	0	0	0	0	2	0
CD8+ T cell	0	12	306	1	0	0	0	0	0	0
gd T cell	0	0	10	37	0	0	0	0	0	0
NKT cell	1	0	0	0	0	6	0	0	0	0
NK cell	0	0	0	0	0	300	0	0	0	0
C Monocyte	1	0	0	0	0	0	1130	14	3	0
NC-Int Monocyte	0	0	0	0	0	0	1	52	1	0
cDC	0	0	0	0	0	0	4	0	76	0
pDC	0	0	0	0	0	0	0	0	3	4", data.table = FALSE
)


precise <- data.table::fread("	B cell	CD4+ T cell	CD8+ T cell	gd T cell	NKT cell	NK cell	C Monocyte	NC-Int Monocyte	cDC	pDC
B cell	332	0	0	0	0	0	5	0	2	0
CD4+ T cell	0	1455	66	3	4	0	0	0	1	0
CD8+ T cell	0	66	277	1	3	4	0	0	0	0
gd T cell	0	5	50	3	0	0	0	0	0	0
NKT cell	0	1	1	0	5	4	0	0	0	0
NK cell	0	0	2	0	3	310	1	0	0	0
C Monocyte	0	0	0	0	0	1	1151	14	4	1
NC-Int Monocyte	0	0	0	0	0	0	3	54	2	0
cDC	0	0	0	0	0	0	17	0	75	0
pDC	1	0	0	0	0	0	0	0	3	18", data.table = FALSE)

or <- c("B cell", "CD4+ T cell", "CD8+ T cell", "gd T cell", 
        "NKT cell", "NK cell", "C Monocyte", "NC-Int Monocyte", "cDC", "pDC")

or_method <- c("harmony", "seurat", "precise")


total <- data.table::fread("339	B cell
1539	CD4+ T cell
363	CD8+ T cell
69	gd T cell
20	NKT cell
337	NK cell
1171	C Monocyte
59	NC-Int Monocyte
92	cDC
22	pDC", data.table = FALSE) %>% set_names(c("total", "V1"))



list(harmony = harmony, seurat = seurat, precise = precise) %>% 
  map(~ pivot_longer(., -V1, names_to = "prediction", values_to = "counts")) %>% 
  bind_rows(.id = "method") %>% 
  filter(V1 == prediction) %>% 
  full_join(total, by = "V1") %>% 
  mutate(acc = round(counts / total * 100, 2)) %>% 
  mutate(V1 = factor(V1, or)) %>% 
  mutate(method = factor(method, or_method)) -> data
  


ggplot(data, aes(V1, acc, colour = method)) +
  geom_linerange(aes(x = V1, ymin = 0, ymax = acc, colour = method), 
                 position = position_dodge(width = 0.5), size = 0.4) +
  geom_point(position = position_dodge(width = 0.5)) +
  xlab("Cell type") +
  ylab("Accuracy") +
  labs(color = guide_legend(title = "Method")) +
  scale_color_viridis_d() +
  theme_classic() +
  rotate_x()



ggplot(data, aes(V1, counts, colour = method)) +
  geom_linerange(aes(x = V1, ymin = 0, ymax = counts, colour = method), 
                 position = position_dodge(width = 0.5), size = 0.2) +
  geom_point(position = position_dodge(width = 0.5)) +
  xlab("Cell type") +
  ylab("Correct cells") +
  labs(color = guide_legend(title = "Method")) +
  scale_color_viridis_d() +
  theme_classic() +
  rotate_x()



