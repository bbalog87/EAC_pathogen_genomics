library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(forcats)

## process Biosample data : extract relevant column

acinetobacter_eac<- acinetobacter_eac_BioSample %>% 
   select(BioSampleAccession, BioSampleBioProjectAccession,
          BioSampleSRAAccession, BioSampleTitle, BioSampleOrganism,
          BioSampleCollectionDate, BioSampleSequencedBy,BioSampleOrganization, 
          BioSampleGeographicLocation)


# Function to extract relevant columns from a table
extract_columns <- function(table) {
  selected_columns <- select(table, BioSampleAccession, 
                             BioSampleBioProjectAccession,
                             BioSampleSRAAccession, 
                             BioSampleTitle, BioSampleOrganism,
                             BioSampleCollectionDate, 
                             BioSampleSequencedBy,
                             BioSampleOrganization, 
                             BioSampleGeographicLocation)
  return(selected_columns)
}

# List of pathogens tables
biosample.list <- list(acinetobacter_eac_BioSample, enterococcus_eac_BioSample,
                       escherichia_eac_BioSample, klebsiella_eac_BioSample,
                       pseudomonas_eac_BioSample, salmonella_eac_BioSample,
                       staphylococcus_eac_BioSample, streptoccocus_eac_BioSample)

# Extract relevant columns and merge all tables
merged_pathognes_eac.df <- bind_rows(lapply(biosample.list, extract_columns))

# Define counbtries list
eac_partnerStates<- c("Burundi", "Kenya", "Tanzania",
                        "Uganda", "Democratic Republic of Congo", "Congo",
                        "DR Congo", "South Sudan","Rwanda")

# Subset the table for EAC partner States
merged_genomes_eac.df <- Bacteria_EAC_genome %>%
  filter(str_detect(Isolation.Country, 
                    regex(paste(eac_partnerStates, collapse = "|"),
                          ignore_case = TRUE))) %>% 
 select(4, 14, 36, 41, 42, 43, 44, 45, 47, Isolation.Country)

bacterial_genomes_eac.df <- BVBRC_Bacteria_Africa %>%
  filter(str_detect(Isolation.Country, 
                    regex(paste(eac_partnerStates, collapse = "|"),
                          ignore_case = TRUE))) %>% 
  select(4, Genus,14,15, Genome.Quality,Completion.Date, BioSample.Accession,
         Assembly.Accession,SRA.Accession, GenBank.Accessions,
         Sequencing.Center,Sequencing.Platform, Sequencing.Depth,
         Contigs, Contig.N50, Coarse.Consistency,Isolation.Source,
         Isolation.Country,Host.Group) %>%
       filter(!is.na(Year)) %>% filter(Year>=2008)



bacterial_genomes_eac.df <- bacterial_genomes_eac.df %>%
  mutate(Pathogens = if_else(
    str_detect(Genus, 
    "Enterococcus|Staphylococcus|Escherichia|Klebsiella|Acinetobacter|Pseudomonas|Salmonella|Streptococcus"),
    "GLASS","Other"))


asm.data <- bacterial_genomes_eac.df %>% 
  filter(Genome.Status != "Plasmid", Assembly.Accession != "") %>% 
  select(Completion.Date, Assembly.Accession) %>%
  mutate(Year = format(as.Date(Completion.Date), "%Y")) %>%
  group_by(Year) %>%
  summarize(AssemblyCount = n()) %>%
  filter(!is.na(Year)) %>%
  filter(Year>=2008)


## Process BVBRC data 
#======================

#### Data from all Africa
all.pathognes <- BVBRC_genome %>%
  select(Assembly.Accession, Genome.Status,
         Completion.Date, Genome.Quality, Geographic.Group) %>%
  filter(Genome.Status %in% c("WGS", "Complete")) %>%
  filter(!is.na(Completion.Date) & Completion.Date != "") %>%
  filter(Geographic.Group!= "") %>%
  mutate(Year = format(as.Date(strptime(Completion.Date, 
                            "%Y-%m-%dT%H:%M:%SZ")), "%Y")) %>%
  filter(Year %in% 2008:2023)

all.pathognes <- all.pathognes %>%
  mutate(Region = ifelse(Geographic.Group == "Oceania" | Geographic.Group == "Antarctica", "Others", Geographic.Group)) %>%
  distinct(Assembly.Accession, .keep_all = TRUE)



summary <- all.pathognes %>%
  filter(Year <= "2016") %>%
  group_by(Region) %>%
  summarize(Total_Assemblies_2016 = n())

summary <- summary %>%
  left_join(all.pathognes %>%
              filter(Year <= "2019") %>%
              group_by(Region) %>%
              summarize(Total_Assemblies_2019 = n()),
            by = "Region") %>%
  left_join(all.pathognes %>%
              filter(Year <= "2022") %>%
              group_by(Region) %>%
              summarize(Total_Assemblies_2022 = n()),
            by = "Region")


assembly_sum <- all.pathognes %>%
  filter(Year >= "2008" & Year <= "2022") %>%
  group_by(Region, Year) %>%
  summarize(Cumulative_Assembly_Count = n()) %>%
  group_by(Region) %>%
  mutate(Cumulative_Assembly_Count = cumsum(Cumulative_Assembly_Count))




assembly_sum$Year <- as.numeric(assembly_sum$Year)  # Convert Year to numeric
# Define the color palette
custom_colors <- c("#FF1F5B",  "#009ADE", "#AF58BA", 
                   "#FFC61E", "#F28522", "#00CD6C")
asm.world<-ggplot(assembly_sum, aes(x = Year,
                         y = Cumulative_Assembly_Count, color = Region)) +
  geom_line(size=2) +
  geom_point(size=4) +
  labs(x = "", y = "Cumulative Assemblies Count", color = "Region",
       title = "Bacterial assemblies globally") +
  scale_x_continuous(breaks = c(2008, 2012,2016,2019,2022)) +
  scale_color_manual(values = custom_colors) +  # Set the custom color palette
  theme_minimal()+
  theme(plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold", color = "black"),
    axis.text = element_text(size = 24, face = "bold", colour = "black"),
    legend.text = element_text(size = 30),
    legend.title = element_text(size=30, face = "bold"),
    legend.key.size = unit(3, "lines"))

ggsave("asembly_count_all3.pdf", width = 13, height = 10, dpi = 6000)

asm.world+eac.asm.plot

### Genomes distribtion EAC

eac.genomes<- BVBRC_genome %>% 
  filter(str_detect(Isolation.Country, 
                   regex(paste(eac_partnerStates, collapse = "|"),
                         ignore_case = TRUE))) %>%
    select(Assembly.Accession, Genome.Status, Genus,
                    Completion.Date, Genome.Quality, 
           Isolation.Country, Host.Group) %>%
  filter(Genome.Status %in% c("WGS", "Complete")) %>%
  filter(!is.na(Completion.Date) & Completion.Date != "") %>%
  filter(Isolation.Country!= "") %>% filter(Genus!= "") %>%
  mutate(Year = format(as.Date(strptime(Completion.Date, 
                                        "%Y-%m-%dT%H:%M:%SZ")), "%Y")) %>%
  filter(Year %in% 2008:2023) %>%
  filter(Isolation.Country!="Congo")


eac.genomes.df <- eac.genomes %>%
  mutate(Pathogens = if_else(
    str_detect(Genus, 
               "Enterococcus|Staphylococcus|Escherichia|Klebsiella|Acinetobacter|Pseudomonas|Salmonella|Streptococcus"),
    "GLASS","Other"))


## Plot donots

color_palette1 <- c("#B7E6A5", "#7CCBA2", "#46AEA0", "#089099", "#003147")



labels <- c("Other", "Avian", "Human", "Nonhuman Mammal", "Plant")
values <- c(2539, 78, 1783, 205, 30)



df <- data.frame(labels, values)
# Calculate total count
total_count <- sum(df$values)

# Calculate percentage
df$percent <- df$values / total_count * 100

# Sort the data frame by the percentage values in descending order
df <- df[order(-df$percent), ]

ggplot(df, aes(x = "", y = values, fill = labels)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_palette) +
  theme_void() +
  theme(legend.position = "none", 
        legend.direction = "horizontal", legend.box = "horizontal",
        legend.margin = margin(10, 0, 0, 0), legend.title = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5, vjust = 1),
        text = element_text(size = 24)) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0), linetype = "dotted", color = "gray") +
  geom_polygon(data = data.frame(x = c(0.3, 0.3, 0.7, 0.7),
                                 y = c(0.3, 0.7, 0.7, 0.3)),
               aes(x, y), color = "black", fill = NA, size = 1) +
  geom_text(aes(y = cumsum(values) - 0.5 * values, label = paste0(round(percent, 1), "%")), size = 5) +
  geom_segment(data = df,
               aes(x = 1.8, y = cumsum(values) - 0.5 * values, xend = 2.5, yend = cumsum(values) - 0.5 * values,
                   color = labels), size = 1.2) +
  geom_text(data = df,
            aes(x = 2.8, y = cumsum(values) - 0.5 * values, label = paste0(labels, " (", round(percent, 1), "%)")),
            hjust = 0, size = 6, color = "black") 







# Data
color_palette2 <- c("#CDE5D2", "#9CCEA7", "red", "#40AD5A", "blue", "#06592A")
labels2 <- c("Burundi", "DR Congo", "Kenya", "Rwanda", "Tanzania", "Uganda")
values2 <- c(2, 54, 1443, 45, 2233, 858)

df2 <- data.frame(labels2, values2)

# Sort the data frame by values in descending order
df2 <- df2 %>% arrange(desc(values2))

# Calculate percentage
df2$percent <- df2$values2 / sum(df2$values2) * 100

# Calculate cumulative percentage
df2$cum_percent <- cumsum(df2$percent)

# Create the donut plot
ggplot(df2, aes(x = 1, y = percent, fill = labels2)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_palette2) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(size = 24, hjust = 0.5, vjust = 1),
        text = element_text(size = 24)) +
  geom_text(aes(y = cum_percent - 0.5 * percent, label = paste0("(", values2, ")")),
            size = 5, color = "black") +
  geom_segment(aes(x = 0.6, y = cum_percent - 0.5 * percent, xend = 0.8, yend = cum_percent - 0.5 * percent),
               size = 1.2, color = "black") +
  geom_text(aes(x = 0.9, y = cum_percent - 0.5 * percent, label = labels2),
            hjust = 0, size = 6, color = "black") +
  xlim(0.5, 1.5) +
  ylim(0, 100) +
  labs(title = "Samples Isolated in East African Countries") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_rect(fill = "white"))




# Data
color_palette2 <- c("#CDE5D2", "#9CCEA7", "red", "#40AD5A", "blue", "#06592A")
labels2 <- c("Burundi",  "Rwanda","DR Congo", "Uganda", "Kenya",  "Tanzania")
values2 <- c(2, 45, 54,858, 1443,2233)

df2 <- data.frame(labels2, values2)

# Sort the data frame by values in descending order
df2 <- df2 %>% arrange(desc(values2))

# Calculate percentage
df2$percent <- df2$values2 / sum(df2$values2) * 100

# Calculate the cumulative percentage
df2$cum_percent <- cumsum(df2$percent) - 0.5 * df2$percent



# Create the donut plot
ggplot(df2, aes(x = 1, y = percent, fill = labels2)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_palette2) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(size = 24, hjust = 0.5, vjust = 1),
        text = element_text(size = 24)) +
  geom_text(aes(y = cum_percent, label = paste0(labels2, " (n = ", values2, ")")), color = "black", size = 5) +
  xlim(0.5, 1.5) +
  ylim(0, 100) +
  labs(title = "Samples Isolated in East African Countries") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_rect(fill = "white"))











library(ggplot2)
library(dplyr)

# Data
labels2 <- c("Burundi", "DR Congo", "Kenya", "Rwanda", "Tanzania", "Uganda")
values2 <- c(2, 54, 1443, 45, 2233, 858)

df2 <- data.frame(labels2, values2)

# Sort the data frame by values in descending order
df2 <- df2 %>% arrange(desc(values2))

# Calculate percentage
df2$percent <- df2$values2 / sum(df2$values2) * 100

# Calculate the cumulative percentage
df2$cum_percent <- cumsum(df2$percent) - 0.5 * df2$percent

# Sequential color palette
color_palette2 <- colorRampPalette(c("#CDE5D2", "#06592A"))(length(labels2))

# Create the donut plot
ggplot(df2, aes(x = 1, y = percent, fill = labels2)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_palette2) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(size = 24, hjust = 0.5, vjust = 1),
        text = element_text(size = 24)) +
  geom_text(aes(y = cum_percent, label = paste0("(", values2, ")")), color = "black", size = 5) +
  xlim(0.5, 1.5) +
  ylim(0, 100) +
  labs(title = "Isolation Country") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_rect(fill = "white"))









df2 <- data.frame(labels2, values2)
df2$percent <- df2$values2 / sum(df2$values2) * 100

df2 <- df2[order(-df2$values2), ]

ggplot(df2, aes(x = "", y = values2, fill = labels2)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_palette2) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(size = 24, hjust = 0.5, vjust = 1),
        text = element_text(size = 24)) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0), linetype = "dotted", color = "gray") +
  geom_polygon(data = data.frame(x = c(0.3, 0.3, 0.7, 0.7),
                                 y = c(0.3, 0.7, 0.7, 0.3)),
               aes(x, y), color = "black", fill = NA, size = 1) +
  geom_text(aes(y = cumsum(values2) - 0.5 * values2, label = paste0("(", values2, ")")), size = 5) +
  #geom_segment(data = df2,
              # aes(x = 1.8, y = cumsum(values2) - 0.5 * values2, xend = 2.5, yend = cumsum(values2) - 0.5 * values2,
                #   color = labels2), size = 1.2) +
  geom_text(data = df2,
            aes(x = 2.8, y = cumsum(values2) - 0.5 * values2, label = paste0(labels2, " (", values2, ")")),
            hjust = 0, size = 6, color = "black") +
  xlim(c(0, 3)) +
  coord_equal(ratio = 1) +
  labs(title = "Samples Isolated in East African Countries") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_rect(fill = "white"))





llibrary(ggplot2)

# Calculate total count
total_count <- sum(df$values)

# Calculate percentage
df$percent <- df$values / total_count * 100

# Create donut plot
plot <- ggplot(df, aes(x = 1, y = percent, fill = labels)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF", "#FFB400", "#A9A9A9")) +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
        legend.margin = margin(10, 0, 0, 0), legend.title = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5, vjust = 1),
        text = element_text(size = 24)) +
  labs(title = "Samples Isolated in the EAC")

# Calculate label positions
df$label_pos <- cumsum(df$percent) - 0.5 * df$percent

# Add label lines and text outside the plot
plot +
  geom_segment(data = df,
               aes(x = 1.5, xend = 1.7, y = label_pos, yend = label_pos),
               color = "black") +
  geom_text(data = df,
            aes(x = 1.8, y = label_pos, label = paste0(labels, " (", round(percent, 1), "%)")),
            hjust = 0, size = 4) +
  geom_circle(x0 = 1, y0 = 0, radius = 0.5, color = "white", fill = "white", size = 4) +
  coord_polar(theta = "y") +
  theme(plot.margin = margin(1, 1, 1, 3, "cm"))  # Adjust the plot margin to accommodate the labels outside









assembly_sum$Year <- as.numeric(assembly_sum$Year)  # Convert Year to numeric

# Define the color palette
custom_colors <- c("#FF1F5B", "#00CD6C", "#009ADE", "#AF58BA", "#FFC61E", "#F28522")

# Custom label function for y-axis
custom_labels <- function(x) {
  labels <- label_number()(x)
  paste0("10^", log10(labels))
}

ggplot(assembly_sum, aes(x = Year,
                         y = Cumulative_Assembly_Count, color = Region)) +
  geom_line() +
  geom_point() +
  labs(x = "Year", y = "Cumulative Assembly Count", color = "Region") +
  scale_x_continuous(breaks = seq(2008, 2022, 2)) +
  scale_y_log10(labels = custom_labels) +  # Add log scale with custom labels
  scale_color_manual(values = custom_colors) +  # Set the custom color palette
  theme_minimal()







# Extract year from Completion.Date column and count assemblies per year
assembly_count <- aggregate(data$Assembly.Accession, by = list(Year = format(as.Date(data$Completion.Date), "%Y")), FUN = length)
colnames(assembly_count) <- c("Year", "AssemblyCount")

# Sample data
data <- data.frame(
  Year = c(2012, 2012, 2014, 2015, 2016),
  AssemblyCount = c(1, 1, 1, 1, 3)
)

# Line plot with marker dot
plot(data$Year, data$AssemblyCount, type = "o", lty = 1, col = "blue", xlab = "Year", ylab = "Assembly Count", main = "Evolution of Assembly Count")
points(data$Year, data$AssemblyCount, col = "blue", pch = 16)

# Add gridlines
grid()

# Customize axis labels
axis(1, at = data$Year, labels = data$Year)
axis(2, at = 0:max(data$AssemblyCount))

# Add legend
legend("topright", legend = "Assembly Count", col = "blue", lty = 1, pch = 16)

# Save the plot as an image (optional)
# dev.copy(png, "assembly_plot.png")
# dev.off()




### Barplots

# Create the dataframe
data <- data.frame(
  Genus = c("Acinetobacter", "Enterococcus", "Escherichia", "Klebsiella", "Pseudomonas", "Salmonella", "Staphylococcus", "Streptococcus", "Non-GLASS"),
  Count = c(138, 28, 891, 232, 100, 241, 220, 441, 2535)
)

# Print the dataframe
data


# Define the color palette
color_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")


library(ggplot2)
library(dplyr)

# Create the dataframe
data <- data.frame(
  Genus = c("Acinetobacter", "Enterococcus", "Escherichia", "Klebsiella", "Pseudomonas", "Salmonella", "Staphylococcus", "Streptococcus", "Non-GLASS"),
  Count = c(138, 28, 891, 232, 100, 241, 220, 441, 2535),
  Type = c(rep("GLASS", 8), "Non-GLASS")
)

# Define the color palette
color_palette <- c("#999999", "#E69F00")  # Colors for GLASS and Non-GLASS





# Create the dataframe
data <- data.frame(
  Genus = c("Acinetobacter", "Enterococcus", "Escherichia", "Klebsiella", "Pseudomonas", "Salmonella", "Staphylococcus", "Streptococcus", "Non-GLASS"),
  Count = c(138, 28, 891, 232, 100, 241, 220, 441, 2535),
  Type = c(rep("GLASS", 8), "Non-GLASS")
)

# Define the color palette
color_palette <- c("#009392", "#E69F00")  # Colors for GLASS and Non-GLASS

# Add line breaks to x-axis labels
data$Genus <- gsub(" ", "\n", data$Genus)

# Create the bar chart

library(ggplot2)
library(dplyr)

# Sort the data by Count in descending order
data <- data %>% arrange(desc(Count))

# Add line breaks to x-axis labels
data$Genus <- gsub(" ", "\n", data$Genus)

# Create the bar chart
ggplot(data, aes(x = reorder(Genus, -Count), y = Count, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_blank(),  # Remove y-axis labels
    axis.title = element_text(size = 14),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray", linetype = "dashed"),
    plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  geom_text(aes(label = paste0("n =", Count)), position = position_stack(vjust = 0.5), size = 12, color = "black") +
  labs(title = "Genome Counts for Each Genus", x = NULL, y = "Count") +
  coord_flip()



library(ggplot2)
library(dplyr)

# Sort the data by Count in descending order
data <- data %>% arrange(desc(Count))

# Add line breaks to x-axis labels
data$Genus <- paste(data$Genus, "\n(n =", data$Count, ")")

# Create the bar chart
ggplot(data, aes(x = reorder(Genus, -Count), y = Count, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_blank(),  # Remove y-axis labels
    axis.title = element_text(size = 24),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray", linetype = "dashed"),
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 1),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 12, color = "black") +
  labs(title = "Assemblies count per GLASS priority pathogens", x = NULL, y = "Count") +
  coord_flip() +
  scale_x_discrete(labels = function(x) gsub("\n", "\n  ", x))  # Add space indentation for x-axis labels



#### Map of seqeuning data EAC

eac.table <- bacterial_genomes_eac.df %>%
  filter(Genome.Status %in% c("WGS", "Complete")) %>%
    select(Completion.Date,Species, Contig.N50,
           Pathogens,Isolation.Country, Genome.Status) %>%
mutate(Year = format(as.Date(strptime(Completion.Date, 
                                      "%Y-%m-%dT%H:%M:%SZ")), "%Y"))
   write.csv(eac.table, "eac.table.txt", row.names = FALSE)
   
  
 eac.table.df <- eac.table %>% filter(Year>2007) %>%
     group_by(Species, Year, Isolation.Country) %>%
     summarise(assembliesCount = n(),
               MaxN50 = max(Contig.N50)) %>%
     ungroup() %>% filter(Species!="")
   
   # Print the resulting dataframe
   print(eac.table.df)
   
   
  
   
   # Your aggregated data
   # Assuming it is stored in a variable called "result"
   
   # Convert Year to numeric for proper ordering on the x-axis
   eac.table.df$Year <- as.numeric(eac.table.df$Year)
   
   # Create the plot
   ggplot(eac.table.df %>% filter(Year>2007),
          aes(x = Year, y = log(MaxN50), size = assembliesCount, color = Pathogens)) +
     geom_point(alpha = 0.7) +
     scale_size_continuous(range = c(3, 20)) +
     labs(title = "Distribution of Assembly Counts over the Years",
          x = "Year", y = "MaxN50") +
     theme_minimal() +
     theme(text = element_text(size = 16),
           legend.position = "right")
   
   
   
   # Create the plot
   ggplot(df, 
          aes(x = Period , y = log(maxN50), size = assembliesCount, color = Pathogens)) +
     geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.7) +
     scale_size_continuous(range = c(3, 20)) +
     labs(title = "Distribution of Assembly Counts over the Years",
          x = "Year", y = "MaxN50") +
     theme_minimal() +
     theme(text = element_text(size = 16),
           legend.position = "right")
   
   


df<-eac.table
   # Create the "Genus" column
   df <- df %>%
     mutate(Genus = sub("^(\\S+).*", "\\1", Species))
   
   # Create the "Period" column
   df <- df %>% filter(Year>2007) %>%
     mutate(Period = case_when(
       Year >= 2008 & Year <= 2012 ~ "2008-2012",
       Year >= 2013 & Year <= 2017 ~ "2013-2017",
       Year >= 2018 & Year <= 2022 ~ "2018-2022"
     )) %>% filter(!is.na(Period)) %>%
     filter(Isolation.Country!="Sudan")
   
   # Group by Genus, Country, and Period and count the number of assemblies and get the maximum N50
 df <- df %>%
     group_by(Genus, Isolation.Country, Period, Pathogens) %>%
     summarize(assembliesCount = n(), maxN50 = max(Contig.N50)) %>%
     ungroup() %>% filter(Genus!="") %>% 
     filter(Genus !="uncultured") %>% 
     filter(Genus !="unicellular")
   
   # Print the updated dataframe
   print(df)
   
   
   
   
   # Load the required library
   library(ggplot2)
   
   # Assuming your data is in a dataframe called "df"
   
   # Create the jitter plot
   ggplot(df, aes(x = Period, y = log10(maxN50), color = Pathogens, shape = Pathogens)) +
     geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.7, size = 4) +
     labs(title = "Distribution of maxN50 by Year and Pathogens Type",
          x = "Year Period",
          y = "maxN50",
          color = "Pathogens Type",
          shape = "Pathogens Type") +
     scale_shape_manual(values = c("GLASS" = 16, "Other" = 17)) +
     theme_minimal() +
     theme(legend.position = "bottom",
           legend.title = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1),
           plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
           text = element_text(size = 12))
   
   
   
   # Create the grouped jitter plot
   ggplot(df, aes(x = Period, y = log10(maxN50), color = Pathogens, shape = Pathogens)) +
     geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.7, size = 4) +
     labs(title = "Distribution of maxN50 by Year and Pathogens Type",
          x = "Year Period",
          y = "log10(maxN50)",
          color = "Pathogens Type",
          shape = "Pathogens Type") +
     scale_shape_manual(values = c("GLASS" = 16, "Other" = 17)) +
     theme_minimal() +
     theme(legend.position = "bottom",
           legend.title = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1),
           plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
           text = element_text(size = 12)) +
     facet_wrap(~ Pathogens)
   
   
   
   
   # Create a new column combining Pathogens and Period
   df$Pathogens_Period <- paste(df$Pathogens, df$Period, sep = " - ")
   
   # Create the grouped jitter plot
   ggplot(df, aes(x = Pathogens_Period, y = log10(maxN50), color = Pathogens, shape = Pathogens)) +
     geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.7, size = 4) +
     labs(title = "Distribution of maxN50 by Year and Pathogens Type",
          x = "Pathogens and Year Period",
          y = "log10(maxN50)",
          color = "Pathogens Type",
          shape = "Pathogens Type") +
     scale_shape_manual(values = c("GLASS" = 16, "Other" = 17)) +
     theme_minimal() +
     theme(legend.position = "bottom",
           legend.title = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1),
           plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
           text = element_text(size = 12)) +
     facet_wrap(~ Pathogens_Period, ncol = 1)
   
   
   # Create the grouped jitter plot
   ggplot(df, aes(x = Period, y = log10(maxN50), color = Pathogens, shape = Pathogens)) +
     geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.7, size = 4) +
     labs(title = "Distribution of maxN50 by Year and Pathogens Type",
          x = "Year Period",
          y = "log10(maxN50)",
          color = "Pathogens Type",
          shape = "Pathogens Type") +
     scale_shape_manual(values = c("GLASS" = 16, "Other" = 17)) +
     theme_minimal() +
     theme(legend.position = "bottom",
           legend.title = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1),
           plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
           text = element_text(size = 12)) +
     facet_grid(Pathogens ~ ., scales = "free_x", space = "free_x", switch = "x")
   

   
   p2 <- df %>%
     drop_na() %>%
     ggplot(aes(x = Period,
                y = log10(maxN50),
                color = Pathogens))+
     geom_boxplot(outlier.shape = NA)+
     geom_point(position = position_jitterdodge(seed = 42))+
     theme(legend.position = "none")+
     scale_color_brewer(palette="Dark2")
   
   
   
   
   
   
   
   df %>%
     drop_na() %>%
     ggplot(aes(x = Period, y = log10(maxN50), color = Pathogens)) +
     geom_boxplot(outlier.shape = NA) +
     geom_point(position = position_jitterdodge(seed = 42)) +
     scale_color_brewer(palette = "Dark2") +
     labs(title = "Distribution of maxN50 by Year and Pathogens Type",
          x = "Year Period",
          y = "log10(maxN50)",
          color = "Pathogens Type") +
     theme_minimal() +
     theme(legend.position = "bottom",
           legend.title = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1),
           plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
           text = element_text(size = 12))

   
   
   
   # Load the required libraries
   
   # Your data (assuming the tibble is named "data")
   # Replace "data" with the actual name of your tibble if it's different
   
   # Create the bubble plot
   ggplot(df, aes(x = Period, y = Isolation.Country, 
                  size = assembliesCount, color = Pathogens)) +
     geom_point(alpha = 0.7) +
     scale_size_continuous(range = c(2, 20)) +
     labs(title = "Bubble Plot of Assemblies by Period, Country, and Pathogens Type",
          x = "Period",
          y = "Isolation Country",
          size = "Assemblies Count",
          color = "Pathogens Type") +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
           text = element_text(size = 12))
   
   
   
   # Create the stacked bar plot
   ggplot(df, aes(x = Isolation.Country, y = assembliesCount, fill = Period)) +
     geom_bar(stat = "identity") +
     facet_grid(. ~ Pathogens) +
     labs(title = "Number of Assemblies per Country and Period (Colored by Pathogens Type)",
          x = "Isolation Country",
          y = "Assemblies Count",
          fill = "Period") +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           strip.text = element_text(size = 12),
           plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
           text = element_text(size = 12))
   
   # Load the required libraries
   library(ggplot2)
   
   # Your data (assuming the tibble is named "data")
   # Replace "data" with the actual name of your tibble if it's different
   
   # Create the grouped bar plot
   ggplot(df, aes(x = Isolation.Country, y = assembliesCount, fill = Pathogens)) +
     geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
     facet_wrap(~ Period, ncol = 3) +
     labs(title = "Number of Assemblies per Country and Period (Colored by Pathogens Type)",
          x = "Isolation Country",
          y = "Assemblies Count",
          fill = "Pathogens Type") +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           strip.text = element_text(size = 12),
           plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
           legend.position = "bottom",
           legend.title = element_blank(),
           text = element_text(size = 12))
   
   
   
   
   
   # Load the required libraries
   library(ggplot2)
   library(dplyr)
   
   # Your data (assuming the tibble is named "df")
   # Replace "df" with the actual name of your tibble if it's different
   
   # Reorder the levels of Isolation.Country based on the number of GLASS pathogens
   df3 <- df %>%
     mutate(Isolation.Country = factor(Isolation.Country, levels = df %>%
                                         group_by(Isolation.Country) %>%
                                         filter(Pathogens == "GLASS") %>%
                                         summarize(count = sum(assembliesCount)) %>%
                                         arrange(count) %>%
                                         pull(Isolation.Country)))
   
   # Create the grouped bar plot with same width bars and sorted by the number of GLASS pathogens
   ggplot(df3, aes(x = Isolation.Country, y = assembliesCount, fill = Pathogens)) +
     geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
     facet_wrap(~ Period, ncol = 3) +
     labs(title = "Number of Assemblies per Country and Period (Colored by Pathogens Type)",
          x = "Isolation Country",
          y = "Assemblies Count",
          fill = "Pathogens Type") +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           strip.text = element_text(size = 12),
           plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
           legend.position = "bottom",
           legend.title = element_blank(),
           text = element_text(size = 12))
   
   
   
   
   
   # Load the required libraries
   library(ggplot2)
   library(dplyr)
   
   # Your data (assuming the tibble is named "df")
   # Replace "df" with the actual name of your tibble if it's different
   
   # Reorder the levels of Isolation.Country based on the number of GLASS pathogens
   df <- df %>%
     mutate(Isolation.Country = factor(Isolation.Country, levels = df %>%
                                         group_by(Isolation.Country) %>%
                                         filter(Pathogens == "GLASS") %>%
                                         summarize(count = sum(assembliesCount)) %>%
                                         arrange(count) %>%
                                         pull(Isolation.Country)))
   
   ggplot(df_summary, aes(x = Isolation.Country,
                          y = Proportion*100, fill = Region)) +
     geom_bar(stat = "identity") +
     labs(title = "Proportion of Pathogen Genomes Sequenced in EAC, Africa, and Out-of-Africa",
          x = "",
          y = "Proportion (%)",
          fill = "Region") +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 30, color = "black"),
           axis.text.y = element_text(size = 30, color = "black"),
           axis.title = element_text(size = 30, color = "black"),
           legend.position = "top",
           legend.title = element_blank(),
           legend.text = element_text(size = 30, color = "black", face = "bold"),
           plot.title = element_text(size = 30, hjust = 0.5, face = "bold"))
   
   
   
   ggplot(df_summary, aes(x = Isolation.Country, y = Proportion * 100, fill = Region)) +
     geom_bar(stat = "identity") +
     labs(title = "Proportion of Pathogen Genomes Sequenced\n in EAC, Africa, and Out-of-Africa",
          x = "",
          y = "Proportion (%)",
          fill = "Region") +
     scale_fill_manual(values = c("#089009", "#E9002D", "#089099")) +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 30, color = "black"),
           axis.text.y = element_text(size = 30, color = "black"),
           axis.title = element_text(size = 30, color = "black"),
           legend.position = "bottom",
           legend.title = element_blank(),
           legend.text = element_text(size = 30, color = "black", face = "bold"),
           plot.title = element_text(size = 30, hjust = 0.5, face = "bold"))
   
   
   
   
   
   
   # Group by Sequencing.Center and calculate the total Number_of_Genomes for each center
   summary_Seq <- world_eac_summary_table %>%
     group_by(Sequencing.Center) %>%
     summarize(Total_Genomes = sum(Number_of_Genomes)) %>%
     arrange(desc(Total_Genomes))
   
   # Print the summary table
   print(summary_Seq)
   
   write.csv(summary_Seq, "summary_Seq.txt", row.names = FALSE)
   
  
   
   
   
   

   
   # Manually set the order of countries
   country_order <- c("Burundi", "Democratic Republic of the Congo", 
                      "Rwanda", "Uganda", "Kenya", "Tanzania")
   
   # Mutate the Isolation.Country column with the desired order
   df <- df %>%
     mutate(Isolation.Country = fct_reorder(Isolation.Country, assembliesCount, .desc = TRUE))
   
   # Create the stacked bar chart
   ggplot(df %>% filter(Isolation.Country != "Congo"), 
          aes(x = Period, y = assembliesCount, fill = Pathogens)) +
     geom_bar(stat = "identity") +
     facet_wrap(~ Isolation.Country, ncol = 3) +
     labs(title = "",
          x = "Period",
          y = "Assemblies Count",
          fill = "Pathogens Type") +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           strip.text = element_text(size = 12),
           plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
           legend.title = element_blank(),
           text = element_text(size = 12))
   
   
   
   

   
   # Manually set the order of countries
   country_order <- c("Rwanda", "Democratic Republic of the Congo", 
                      "Burundi", "Uganda", "Kenya", "Tanzania")
   
   # Mutate the Isolation.Country column with the desired order
   df <- df %>%
     mutate(Isolation.Country = fct_relevel(Isolation.Country, country_order))
   
   # Create the stacked bar chart
   ggplot(df %>% filter(Isolation.Country != "Congo"), 
          aes(x = Period, y = assembliesCount, fill = Pathogens)) +
     geom_bar(stat = "identity") +
     facet_wrap(~ Isolation.Country, ncol = 3) +
     labs(title = "",
          x = "",
          y = "Assemblies Count",
          fill = "Pathogens Type") +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           strip.text = element_text(size = 12),
           plot.title = element_text(size = 16, hjust = 0.5, vjust = 1),
           legend.title = element_blank(),
           text = element_text(size = 12))+
     theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 25, color = "black", face = "bold"),
           strip.text = element_text(size = 30, color = "black", face = "bold"),
           plot.title = element_text(size = 30, hjust = 0.5, vjust = 1, color = "black", face = "bold"),
           legend.title = element_blank(),
           legend.text = element_text(size = 30, color = "black", face = "bold"),
           text = element_text(size = 30, color = "black", face = "bold"))
   
   
   
   
   ggsave("Glass_EAC.pdf", width = 10, height = 10, dpi = 6000)
   
   
   
      

   
   # Filter by sequencing platform
   world.eac$Technology <- case_when(
     is.na(world.eac$Sequencing.Platform) | world.eac$Sequencing.Platform == "" ~ "no.info",
     grepl("PacBio|Sequel", world.eac$Sequencing.Platform, ignore.case = TRUE) ~ "PacBio",
     grepl("Oxford|ONT|Nanopore|ION", world.eac$Sequencing.Platform, ignore.case = TRUE) ~ "ONT",
     grepl("Illumina", world.eac$Sequencing.Platform, ignore.case = TRUE) ~ "Illumina",
     grepl("PacBio", world.eac$Sequencing.Platform, ignore.case = TRUE) & grepl("Illumina|Nanopore", world.eac$Sequencing.Platform, ignore.case = TRUE) ~ "mixed",
     TRUE ~ "Other"
   )
   
   # Print the updated data frame with the 'Technology' column
   print(world.eac)
   write.csv(world.eac, "world.eac.txt", row.names = FALSE)
   
   world.eac_summary_table <- world.eac %>%
     group_by(Sequencing.Center, Isolation.Country, Pathogens) %>%
     summarise(Number_of_Genomes = n())
   
   # Print the summarized table
   print( world.eac_summary_table)
   write.csv(world.eac_summary_table, "world.eac_summary_table.txt",
             row.names = FALSE)
   
      

   result.world.eac <- world_eac_summary_table %>%
     group_by(Owner.Country, ISO3, Pathogens) %>%
     summarize(
       Total_Genomes = sum(Number_of_Genomes),
       GLASS_Count = sum(Number_of_Genomes[Pathogens == "GLASS"]),
       Other_Count = sum(Number_of_Genomes[Pathogens == "Other"])
     )
   


   

   
   # Group by Isolation.Country and Region, and summarize the data
   df_summary <-world_eac_summary_table %>%
     group_by(Isolation.Country, Region) %>%
     summarise(Total_Genomes = sum(Number_of_Genomes)) %>%
     group_by(Isolation.Country) %>%
     mutate(Proportion = Total_Genomes / sum(Total_Genomes)) %>%
     ungroup()
   
   # Create the stacked bar chart
   ggplot(df_summary, aes(x = Isolation.Country, y = Proportion*100, fill = Region)) +
     geom_bar(stat = "identity") +
     labs(title = "Proportion of Pathogen Genomes Sequenced in EAC, Africa, and Out-of-Africa",
          x = "",
          y = "Proportion (%)",
          fill = "Region") +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 30, color = "black"),
           axis.text.y = element_text(size = 30, color = "black"),
           axis.title = element_text(size = 30, color = "black"),
           legend.position = "top",
           legend.title = element_blank(),
           legend.text = element_text(size = 30, color = "black", face = "bold"),
           plot.title = element_text(size = 30, hjust = 0.5, face = "bold"))





# Calculate the total number of genomes for each country and region
total_genomes <- world_eac_summary_table %>%
  group_by(Isolation.Country, Region) %>%
  summarize(Total_Genomes = sum(Number_of_Genomes))

# Calculate the total number of genomes for each country
country_total_genomes <- total_genomes %>%
  group_by(Isolation.Country) %>%
  summarize(Total_Country_Genomes = sum(Total_Genomes))

# Merge the data to get the final summary
summary_data <- merge(total_genomes, country_total_genomes,
                      by = "Isolation.Country")

# Calculate percentages for each category in each country
summary_data <- summary_data %>%
  mutate(
    Percentage_EAC = ifelse(Region == "EAC", Total_Genomes / Total_Country_Genomes * 100, 0),
    Percentage_Africa = ifelse(Region == "Africa", Total_Genomes / Total_Country_Genomes * 100, 0),
    Percentage_Out_of_Africa = ifelse(Region == "Out-of-Africa", Total_Genomes / Total_Country_Genomes * 100, 0)
  ) %>%
  select(-Total_Genomes, -Total_Country_Genomes)

# Display the final summary
print(summary_data)

write.csv(summary_data, "summary_data_sequencing.txt", row.names = FALSE)




# Group the data by Isolation.Country and calculate the sum of Number_of_Genomes
eac_summary <- world_eac_summary_table %>%
  group_by(Isolation.Country, Pathogens ) %>%
  summarize(Total_Genomes = sum(Number_of_Genomes))

# Calculate the percentage of GLASS pathogens for each country
summary_data <- summary_data %>%
  mutate(Percentage_GLASS = (Total_Genomes / sum(Total_Genomes)) * 100)

