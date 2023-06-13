from collections import namedtuple

debugging=False
# %%

JsonComparison = namedtuple("JsonComparison", ["is_identical", "keys", "reason", "left", "right"])

def compare_jsons(
        left, right, stop_on_first_diff=False, compare_func=None, 
        raise_exception_for_unsupported_type = True, keys=[]
    ):
    """Compares two supposed JSONs made exclusively of types [bool, int, float, str, list, dict]

    Does a depth first search to find non-matching values.
    Note: as Sequences, only lists are supported (no tuple, set, range, ...)

    arguments:
    - left, right: jsons to be compared
    - stop_on_first_diff: whether to stop once a single difference has been found (=whether to simply check for any difference or find all differences)
    - compare_func (optional): a custom comparison function taking two arguments and returning a boolean
    - raise_exception_for_unsupported_type: whether to raise an exception if any entry is not of type [bool, int, float, str, list, dict]
    - keys: do not set, used for recursion, indicates which keys have been traversed
    
    returns: a list of JsonComparison namedtuple, each with elements:
    - is_identical: whether the two JSONs are identical
    - keys: the list of keys indicating where the equality test first failed.
    - reason: why the test failed
    - left: the left value where the test failed.
    - right: the right value where the test failed.
    If the two JSONs are identical, the list will consist of a single JsonComparison with is_identical=True.
    If they are different and stop_on_first_diff==True, the list will consist of a single JsonComparison for the first place where the test failed.
    If they are different and stop_on_first_diff==False, it will return one JsonComparison per difference.
    """
    if debugging:
        print("compare_jsons() key: ", keys, "left:",left, "right:",right, "left==right:",left==right)
    # simple equality test for int, float, str and bool
    if isinstance(left, float) or isinstance(left, int) or isinstance(left, str) or isinstance(left, bool) or left is None:
        if compare_func is None:
            result = [JsonComparison(left==right, keys, "equality test", left, right)]
            return result
        else:
            return [JsonComparison(compare_func(left, right), keys, "compare_func test", left, right)]
    # if not simple type, ensure types are identical
    elif type(left) != type(right):
        return [JsonComparison(False, keys, "different types", left, right)]
    # dict comparison
    elif isinstance(left, dict):
        if sorted(left.keys()) != sorted(right.keys()):
            return [JsonComparison(False, keys, "dicts with different keys", left, right)]
        results = []
        both_identical = True
        for i,k in enumerate(left.keys()):
            new_results = compare_jsons(left[k], right[k], stop_on_first_diff, compare_func, raise_exception_for_unsupported_type, keys+[k])
            both_identical = both_identical and new_results[0].is_identical
            if stop_on_first_diff and not both_identical:
                return new_results
            results+=new_results
        if both_identical:
            return [JsonComparison(True, keys, "dict equality success", left, right)]
        else:
            return [r for r in results if not r.is_identical]
    # list comparison
    elif isinstance(left, list):
        if len(left) != len(right):
            return [False, keys, "lists with different len()", left, right]
        results = []
        both_identical = True
        for i, lv in enumerate(left):
            new_results = compare_jsons(lv, right[i], stop_on_first_diff, compare_func, raise_exception_for_unsupported_type, keys+[i])
            both_identical = both_identical and new_results[0].is_identical
            if stop_on_first_diff and not both_identical:
                return new_results
            results+=new_results
        if both_identical:
            return [JsonComparison(True, keys, "list equality success", left, right)]
        else:
            return [r for r in results if not r.is_identical]
    # other types are unsupported -> error or False
    else:
        if raise_exception_for_unsupported_type:
            raise Exception(
                f"compare_jsons(): type(left)=={type(left)} not in [bool, int, float, str, list, dict, NoneType]."+ 
                f"\n\ttype(left)={type(left)}, left={left}"+ 
                f"\n\ttype(right)={type(right)}, right={right}"
            )
        else:
            return [JsonComparison(False, keys, "unsupported types", left, right)]


def compare_jsons_types(left, right, **kwargs):
    """compare_json() to compare for types
    
    compare_func = lambda x,y: type(x)==type(y)
    note: unlike JSON, float and int types are distinct
    """
    return compare_jsons(left, right, compare_func=lambda x,y: type(x)==type(y), **kwargs)

def compare_jsons_types_json_numbers(left, right, **kwargs):
    """compare_jsons_types() considering float and ints identical, as in proper JSON"""
    def compare_func(left, right):
        if isinstance(left, float) or isinstance(left, int):
            if isinstance(right, float) or isinstance(right, int):
                return True
        return type(left)==type(right)
    return compare_jsons(left, right, compare_func=compare_func, **kwargs)

def compare_jsons_types_inty_floats(left, right, **kwargs):
    """compare_jsons_types() with handling the case for integer-like floats"""
    def compare_func(left, right):
        left_types = (int, float) if isinstance(left, float) and left == int(left) else (type(left),)
        right_types = (int, float) if isinstance(right, float) and right == int(right) else (type(right),)
        return any(lt in right_types for lt in left_types)
    return compare_jsons(left, right, compare_func=compare_func, **kwargs)

# %%

if False:
    a= {
        "a": True, "b":42, "c":[4,"x",5,{"x":"xy"}], "d":{
            "bb": "bb22", "aa": 1/3
        }
    }
    b= {
        "b":42, "a": True, "c":[4,"x",5,{"z":"xy"}], "d":{
            "aa": 10/3132, "bb": "bb22"
        }
    }

    compare_jsons(a, b)
# %%
