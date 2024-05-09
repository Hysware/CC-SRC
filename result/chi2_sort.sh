if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <column_number>"
    exit 1
fi

input_file="$1"
column_number="$2"

# 检查文件是否存在
if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found."
    exit 1
fi

# 临时文件用于存储奇数行和偶数行的配对数据
temp_file=$(mktemp)

# 提取奇数行和偶数行的数据并配对存储
awk -v col="$column_number" '{
    if (NR % 2 == 0) {
        even_row = $col
        print prev_odd_row "\t" even_row
    } else {
        prev_odd_row = $col
    }
}
END {
    # 如果文件行数为奇数，输出最后一行的奇数行数据
    if (NR % 2 != 0) {
        print prev_odd_row
    }
}' "$input_file" > "$temp_file"

# 按照偶数行的数据进行排序，并输出排序后的偶数行和对应的奇数行
sort -k2 -o "$temp_file" "$temp_file"
awk '{print $1 "\t" $2}' "$temp_file"

# 删除临时文件
# rm "$temp_file"

