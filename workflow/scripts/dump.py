# file to dump code

##################################
##################################
##################################

# retrieve affected genes from gene sets

# metric = "logFC"
# lib = "KEGG"
# lib = "GO"
# tool = "gseapy"
# tool = "string"
# #tool = "clusterProfiler"

# gene_table_file = f"../../results/{project_name}/gene_converter.csv"
# converter = pd.read_csv(gene_table_file, index_col=0)

# terms_dict = dict()

# #terms = summary_dict[lib]["depth_df"].iloc[0]["Description"]

# depth_df = pd.read_csv(f"{savepath}/combined/syn.depth.{lib}.{project_name}.csv", index_col=0)
# d = depth_df
# d = d[(d["Combined FDR"]<0.05) & (d["Depth"]>1)]
# terms = d["Description"].values.tolist()

# #terms = "E-box binding"

# if not isinstance(terms, list):
#     terms = [terms]

# terms = [t.strip() for t in terms]
# terms.append("NADH dehydrogenase activity")

# for t in terms:
#     terms_dict[t] = None

# tab = pd.read_csv(f"{savepath}/syn.{tool}.{metric}.{lib}.{project_name}.csv", index_col=0)

# not_found = []
# for t in terms_dict:
#     print(t)
#     match tool:
#         case "gseapy":
#             lead_genes = tab.loc[t, "Lead_genes"]
#             lead_genes = lead_genes.split(";")
#         case "clusterProfiler":
#             tt = tab[tab["Description"].str.lower()==t.lower()].index
#             lead_genes = tab.loc[tt, "core_enrichment"]
#             lead_genes = lead_genes.values[0].split("/")

#         case "string":
#             tt = tab[tab["Description"].str.lower()==t.lower()].index

#             if "matching proteins in your input (labels)" in tab.columns:
#                 lead_genes = tab.loc[tt, "matching proteins in your input (labels)"]
#             else:
#                 lead_genes = tab.loc[tt, "matching proteins in your input (IDs)"]

#             try:
#                 lead_genes = lead_genes.values[0].split(",")
#             except IndexError:
#                 print("Not found:", t)
#                 not_found.append(t)
#                 #display(tab.head())
#         case default:
#             assert 0
#     terms_dict[t] = lead_genes

# for n in not_found: del terms_dict[n]

# dea = pd.read_csv("../../resources/Chiara/edger.qlf.lfc0.SA_WT.p1.csv", index_col=0)


# ##############

# dd = depth_df[depth_df["Description"].str.contains("dehydrogenase")]
# for d in dd.index:
#     print(depth_df.loc[d,"Description"], depth_df.loc[d,"Configurations"])

# dd

# ##############



# from utils import ensp_to_gene_symbol
# x = pd.read_csv("../../results/Chiara.QLF.SA_WT/proteins/chiara_list.csv")
# x.columns = [x.strip().lower().replace(" - ","-") for x in x.columns]


# for t in terms_dict:
#     lead = terms_dict[t]

#     if "ENSP0" in lead[0]:
#         species = lead[0].split(".")[0]
#         print(species)
#         ensp = [le.split(species)[1][1:] for le in lead]
#         converted = ensp_to_gene_symbol(ensp, species).values
#     elif tool == "clusterProfiler":
#         converted = converter[converter["ENTREZID"].astype(str).isin(terms_dict[t])]["SYMBOL"].values
#     else:
#         converted = lead

#     d = dea[dea["gene_name"].isin(converted)]
#     a = set(dea[dea["gene_name"].isin(converted)]["gene_name"])
#     b = set(terms_dict[t])
#     c = set(x[t.lower()].dropna())

#     if "dehydrogenase" in t:
#         # print(t, f"\nGenes in enrichment table: {len(terms_dict[t])}", f"\nedgeR table: {len(d)}")
#         # print("Difference edger vs enrichment:", a.symmetric_difference(b))
#         print(t, f"\nGenes in enrichment table: {len(terms_dict[t])}")
#         print(f"Genes in Chiara's excel file: {len(c)}")
#         print("Intersection enrichment & excel:", len(b.intersection(c)))
#         print("Union enrichment & excel:", len(b.union(c)))
#         print("In enrichment, missing in excel:", b.difference(c))
#         print("In excel, missing in enrichment:", c.difference(b))
#         print("\n============")

#     d.to_csv(f"{savepath}/proteins/{t}.csv")

##################################
##################################
##################################