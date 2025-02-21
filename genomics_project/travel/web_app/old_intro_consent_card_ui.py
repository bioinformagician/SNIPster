                ui.card(
                    {"class": "centered-container"},
                    
                    ui.h3("Introduction & Consent", class_="card-title"),
                    ui.div(
                        {"class": "scrollable-text"},  # Add this class to make it scrollable
                        ui.div(
                            {"class": "info-text"},
                            ui.tags.head(
                                ui.tags.script(src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"),
                                ui.tags.link(
        rel="stylesheet",
        href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"
    ),
                                
                            ),
                            ui.p(
                                """Our platform provides genetic analysis using polygenic risk score(PGS) to help you understand your health predispositions 
                                and genetic traits. Using advanced genomic analysis techniques, we analyze your unique genetic data to provide 
                                insights about your:"""
                            ),
                            ui.tags.ul(
                                {"class": "feature-list"},
                                ui.tags.li("Disease risk assessments and health predispositions"),
                                ui.tags.li("Inherited traits and characteristics"),
                                ui.tags.li("Athletic performance potential"),
                                ui.tags.li("Uniqueness compared to other users"),
                            ),
                            ui.p(
                                """Simply upload your genetic data file, and our platform will analyze various genetic markers to 
                                provide you with detailed insights about your genetic profile. All analysis is performed with strict 
                                privacy measures in place, meaning that your data cannot be traced back to you. Your genetic data
                                holds a lot of information about you, but in reality the amount of your genome contained in your 
                                file from 23andme, ancestry, etc. Is actually just around 0.0002% of your genome - meaning information about 
                                600.000 base pairs out of the 3.100.000.000 it takes to create exactly you."""
                            ),
                            ui.p(
                                """SNPSTER's PGS approach is to an extend inspired by the previous open source service impute.me which
                                you can learn more about in detail following this PMID 32714365 on pubmed. The formula used for the PGS
                                calculation is given as such:"""
                            ),
                            ui.p([
                                "Population SNP Score: ",
                                ui.tags.span("$$\\text{Population score}_{snp} = \\text{frequency}_{snp} \\times 2 \\times \\text{beta}_{snp}$$"),
                                "Zero-centered Score: ",
                                ui.tags.span("$$\\text{Zero-centered score} = \\sum{\\text{Beta}_{snp} \\times \\text{Effect-allele-count}_{snp} - \\text{Population score}_{snp}}$$"),
                                "Z-score: ",
                                ui.tags.span("$$\\text{Z-score} = \\frac{\\text{Zero-centered score}}{\\text{Standard deviation}_{population}}$$")
                            ]),
                            ui.p(
                                """The data used for the PGS calculation is mostly pulled from the gwas catalog at https://www.ebi.ac.uk/gwas/ 
                                the 20th of november 2024 with some minor contribution from data on SNPEDIA. The 1000 genomes project phase 3 
                                30x dataset (GRCh38) is also being used to impute missing genotypes from raw files"""
                            )
                        )
                    ),
                    ui.div(
                        {"class": "consent-text"},
                        ui.h4("Privacy and Consent"),
                        ui.p(
                            "By uploading your data, you consent to our analysis of your genetic information. ",
                            "Your data will be processed securely and confidentially with no personally identifiable information attached to your information.",
                            " For the purpose of further developing this platform's accuracy and features, your raw file input will be saved in a database with",
                            " no attachment to you identity. This also means that it will be impossible to remove your data on request, as it cannot be traced",
                            " back to you. For any questions, feedback or requests, please contact snpster@gmail.com"
                        ),
                        ui.input_checkbox("consent", "I agree to the privacy, use, and consent terms",  width='100%'),
                    ),
                    ui.h3("Upload Your Genetic Data", class_="card-title"),
                    ui.div(
                            ui.p("Upload your genetic data file to begin the analysis", class_="upload-instruction"),
                            ui.div({"class": "required-field"}, "* Required for analysis")
                        ),
                    ui.input_file(
                        "user_file",
                        "",
                        button_label="Browse Files",
                        placeholder="Drag and drop your file here",
                        multiple=False,
                        accept=".txt,.csv,.tsv",
                        width="100%"
                    ),
                    
                        ui.h4("Select your superpopulation", class_="card-title"),
                        ui.input_radio_buttons(
                            "user_superpopulation",
                            "",  
                            choices=["European", "East Asian", "African", "South Asian", "Admixed American"], 
                            selected="European",
                            width="100%",
                            inline=True
                        ),

                    ui.input_task_button("run_analysis", "Analyze My Genome!", class_="btn btn-primary", width="100%"),
                ),