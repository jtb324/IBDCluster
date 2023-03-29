from dataclasses import dataclass
from typing import Protocol


class FileIndices(Protocol):
    def get_haplotype_id(IID: str, hap_id: str) -> str:
        ...


@dataclass
class HapIBD(FileIndices):
    id1_indx: int = 0
    hap1_indx: int = 1
    id2_indx: int = 2
    hap2_indx: int = 3
    chr_indx: int = 4
    str_indx: int = 5
    end_indx: int = 6
    cM_indx: int = 7

    def get_haplotype_id(IID: str, hap_id: str) -> str:
        return "{0}.{1}".format(IID, hap_id)


@dataclass
class Germline(FileIndices):
    id1_indx: int = 0
    hap1_indx: int = 1
    id2_indx: int = 2
    hap2_indx: int = 3
    chr_indx: int = 4
    str_indx: int = 5
    end_indx: int = 6
    cM_indx: int = 10
    unit: int = 11

    def get_haplotype_id(IID: str, hap_id: str) -> str:
        return hap_id


@dataclass
class iLASH(FileIndices):
    id1_indx: int = 0
    hap1_indx: int = 1
    id2_indx: int = 2
    hap2_indx: int = 3
    chr_indx: int = 4
    str_indx: int = 5
    end_indx: int = 6
    cM_indx: int = 10

    def get_haplotype_id(IID: str, hap_id: str) -> str:
        return hap_id


@dataclass
class Rapid(FileIndices):
    id1_indx: int = 1
    hap1_indx: int = 3
    id2_indx: int = 2
    hap2_indx: int = 4
    chr_indx: int = 0
    cM_indx: int = 7
    str_indx: int = 5
    end_indx: int = 6

    def getHAPID(IID, hapID) -> str:
        return "{0}.{1}".format(IID, hapID)


def create_indices(ibd_file_format: str) -> FileIndices:
    """Factory method to generate the proper file indice object based on the ibd program

    Parameters
    ----------
    ibd_file_format: str
        string indicating what ibd program was used identify IBD segments. EX: hapibd,
        ilash, rapid, and germline. expects this value to be lower case

    Returns
    -------
    FileIndices
        returns an object that conforms to the FileIndices protocol. It will have the
        method getHAPID. It will also have the correct indices for the ibd program

    Raises
    ------
    ValueError
        Raises a value error if the user passes an ibd_file_format that is not hapibd,
        hap-ibd, germline, ilash, rapid
    """
    format_selector = {
        "germline": Germline(),
        "ilash": iLASH(),
        "hapibd": HapIBD(),
        "rapid": Rapid(),
    }

    return format_selector.get(ibd_file_format)
