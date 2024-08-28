# Conservation policies for Marine Megafauna

## Data

MegaMove Data: Is a collection of position registers of tagged individuals of different individuals from Marine Megafauna species. Each entry of the data is a different record, indicationg the individual ID, its species, the time of recording and the GPS location (lat-lon coordinates).


Following the data, individual trajectories are obtained. For each individual, a trajectory is a sequence of points separated by less than 2 days. Each individual may describe multiple trajectories. All trajectories are scaled to last 30 days, so the scale indicates _"how many days per month did the individual spent at some point"_.


### Characterize EEZs
EEZs (Exclusive economic zones) represent areas of eexclusive use of a sovereign country/international entity.
EEZs are coded by an integer in alphabetical order, starting with 0 for High Seas (referred as "-1").

World is divided into $0.5"\times0.5"$ cells, enumerated from the bottom left square. Each cell is designated to the EEZ that occupy the largest fraction of the cell. Geographic records of the data are assigned to an EEZ using these cells:

<center>lat-lon coords. &rarr; cell id &rarr; EEZ</center>

There are economic data of 218 countries according to the **world bank**, and 261 EEZ zones. There are some discrepancies between these datasets:

- 7 countries have no EEZ
- 15 EEZ Zones are conflict zones between different countries.
- The economies of Macao (MAC) and Hong Kong (HKG) do not have EEZs, instead they are within the EEZ of China (CHN).
- Kosovo (XKX) does not have an EEZ, but is within the EEZ of Serbia (SRB).

In total 210 EEZ map directly to a economy. 
- The EEZ of MNP++ is divided in the economies of MNP and GUM, that are both high income economies. *manualy introduced*
- The EEZs of Guemsey and Jersey both belong to the Channel islands economy, that is high income. *manualy introduced*
- Curaçao has two different codes on each dataset. *manualy introduced*
This adds 4 more EEZ to income data, total of 214. 15 EEZs are conflict zones that are not mapped to any economy.

The other 32 EEZs are not mapped to any economy, but are not conflict zones. They are:
- Anguilla (AIA): UK overseas territory
- Netherlands Antilles (ANT): The EEZ contains the island of Sant Eustatius (Politically part of the BES Islands), that is part of the Netherlands (NLD)
- Antartica (ATA): Antartica
- French Southern & Antarctic Lands (ATF): French overseas territory
- Bonaire, Sint Eustatius and Saba (BES): Dutch overseas territory. The EEZ only contains Bonaire and Saba.
- Bouvet Island (BVT): Norwegian overseas territory. Inhabited.
- Cocos (Keeling) Islands (CCK): Australian overseas territory. 
- Cook Islands (COK): New Zealand associated state.
- Clipperton Island (CPT): French overseas territory.       
- Christmas Island (CXR): Australian external territory (non-self-governing).
- Western Sahara (ESH): Sahrawi Arab Democratic Republic territory occupied by morocco.
- Falkland Islands (Malvinas) (FLK): British overseas territory.
- Guadalupe (GLP): French overseas territory.
- French Guiana (GUF): French overseas territory.
- Heard Island and McDonald Islands (HMD): Australian external territory. Inhabited.    
- British Indian Ocean Territory (IOT): British overseas territory. Military base.
- Montserrat (MSR): British overseas territory.
- Martinique (MTQ): French overseas territory.
- Mayotte (MYT): French overseas territory.
- Norfolk Island (NFK): Australian external territory (non-self-governing).
- Niue (NIU): New Zealand associated state.
- Pitcairn (PCN): British overseas territory.
- Reunion (REU): French overseas territory.
- South Georgia and the South Sandwich Islands (SGS): British overseas territory.
- Saint Helena, Ascension and Tristan da Cunha (SHN): British overseas territory.   
- Svalbard (SJM): Norwegian overseas territory.
- Jan Mayen (SJM): Norwegian overseas territory. No permanet populaition.
- Saint Pierre and Miquelon (SPM): French overseas territory.
- Tokelau (TKL): New Zealand dependent territory.
- United States Minor Outlying Islands (UMI): US overseas territory.
- Vatican City (VAT): Vatican City State. 
- Wallis and Futuna (WLF): French overseas territory.


World Bank data splits countries in income quartiles.

After mapping all trajectory points of the dataset to an EEZ, we compute the average timestay of each EEZ. The timestay is defined as the time the individual spent around the registered point, and is computed as half of the time spent from/to the last/next point (if any).  

1. For each individual, aggregate all trajectory points within the same EEZ by summing their time stays. 

2. For each species, sum the timestays of individuals at every EEZ and divide by the number of registered individuals of that species. 

3. For each EEZ, sum the average time stay of all species and divide by the number of species that visited that zone.

As a result, for each EEZ we have computed the **Mean average time stay of the species that visit it**.

***In the first plot***, we show all this information by coloring the marine areas of the EEZ by the average species timestay, the number of species that appear at each EEZ by a circle at the center of the continental part and the economic group by a color of the continental part of the EEZ.

### Network construction

Between individuals and EEZ we define a bipartite network, with links pointing EEZs visited by each individual. In a weighted version, each link is weighted by the timestay of the individual and the EEZ.

The individual-EEZ bipartite network can be aggregated into a species-EEZ network.

Bipartite network can be projected into a EEZ-EEZ network, where links between EEZs account for the number of individuals that visit both zones. A community detection algorithm (Infomap) is run in this network to find communities of EEZs. Results are different if we inlcude High Seas links or not.

_communities could be further investigated, showing the distribution of income group, number of species/individuals, distribution of link propertries (degree, weights) of its components_

### Percolation

There are several ways to compute how to protect the countries

1. Ascending/Descending order of the number of individuals. With/without rich protecting from the beggining.

2. Protecting first those areas that contains the individuals that are easier to protect.


### Nestedness




