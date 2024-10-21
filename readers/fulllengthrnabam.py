import pandas as pd
import numpy as np
import pysam
from bugv.readers.rnabam import _RNABamOnlyReads, convert_reads_to_features
from bugv.readers.feature import Feature, FeatureArray

class RNABamFeature(_RNABamOnlyReads):
        
    def _pre_process_feature(self, feature_array):
        pass
        
    def _load_region_data_func(self):
        reads = self.bam_obj.fetch(self.chr_name,
                                   self.region_start - 1,
                                   self.region_end)
        region_type_name = "mRNA"
        feature_array = convert_reads_to_features(reads, region_type_name, self._get_strand_func)
        self._pre_process_feature(feature_array)
        return feature_array

class FulllengthBam(RNABamFeature):
    
    """
    The tags extract from bam are stored in feature_array.df, but not in each read feature.
    """
    
    DEFAULT_CONFIG = {"min_polya_length": 15}
    
    def __init__(self, file_name="", add_polyA_feature_flag=True, min_polya_length=-1):
        self.file_name = file_name
        self.add_polyA_feature_flag = add_polyA_feature_flag
        if min_polya_length == -1:
            self.min_polya_length = self.DEFAULT_CONFIG["min_polya_length"]
    
    def _get_strand_func(self, read):
        read_strand = ""
        try:
            read_strand = read.get_tag("rs")
        except:
            pass
        if not read_strand:
            read_strand = "-" if read.is_reverse else "+"
        return read_strand
    
    def _pre_process_feature(self, feature_array):
        
        def get_and_set_tags(feature_array,
                             keys=["mi", "ty", "rs", "ir", "pa", "pn"],
                             column_names=["mRNA_id", "read_type", "rna_strand", "intron_retention", "polya_length", "pos_not_consistent"],
                             default_values=["", "unknown", "", False, 0, 1]):
            
            def _get_tag(read, key, default_value):
                # Try to get the tag from the BAM file, otherwise use the default value
                try:
                    if key == "pa":
                        # Adapted to use the 'pt' tag for polyA tail length from Dorado basecaller
                        return read.get_tag("pt")
                    return read.get_tag(key)
                except KeyError:
                    # Special handling for certain tags if they are missing
                    if key == "mi":
                        # Infer mi from the reference name if possible
                        return read.reference_name.split("|")[0] if "|" in read.reference_name else ""
                    elif key == "ir":
                        # Infer intron retention based on the gene type in the reference name
                        gene_type = read.reference_name.split("|")[-1] if "|" in read.reference_name else ""
                        return "intron_retention" in gene_type
                    elif key == "rs":
                        # Use strand information if available
                        return "-" if read.is_reverse else "+"
                    elif key == "pa":
                        # Default polyA length if not present
                        return 0
                    else:
                        return default_value

            list_of_column_data = [[] for _ in range(len(column_names))]
            
            for feature in feature_array.features:
                read = feature.parm["read"]
                for key, default_value, column_name, column_data in zip(keys, default_values, column_names, list_of_column_data):
                    value = _get_tag(read, key, default_value)
                    column_data.append(value)
                    setattr(feature, column_name, value)
                    
            for column_name, column_data in zip(column_names, list_of_column_data):
                feature_array.df[column_name] = column_data

        def add_polyA_feature(region_data, min_polya_length=15):
            for feature in region_data:
                if feature.polya_length >= min_polya_length:
                    if feature.strand == "+":
                        feature.children.append(Feature(feature.end + 1, 
                                                        feature.end + feature.polya_length, 
                                                        feature.chr_name,
                                                        feature.strand, "", 1, "polyA"))
                        feature.end = feature.end + feature.polya_length                 
                    else:
                        feature.children.insert(0, Feature(feature.start - feature.polya_length, 
                                                           feature.start - 1, 
                                                           feature.chr_name,
                                                           feature.strand, "", 1, "polyA"))
                        feature.start = feature.start - feature.polya_length            
            return region_data

        # Process tags and add polyA feature if necessary
        get_and_set_tags(feature_array)
        if self.add_polyA_feature_flag:
            add_polyA_feature(feature_array, self.min_polya_length)
