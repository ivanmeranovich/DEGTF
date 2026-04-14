import sys
import pandas as pd
import os

def clean_fdr(value):
    """
    Вспомогательная функция для очистки значений FDR.
    Преобразует строки с запятыми (например, '0,05') в числа с плавающей точкой.
    """
    if pd.isna(value):
        return None
    try:
        if isinstance(value, str):
            value = value.replace(',', '.')
        return float(value)
    except ValueError:
        return None

def build_deg_mapping(df_deg):
    """
    Создает словарь для быстрого поиска совпадений по таблицам DEG.
    Ключ: target_id или locus (в виде строки)
    Значение: кортеж (logFC, FDR)
    """
    mapping = {}
    for index, row in df_deg.iterrows():
        fdr_val = clean_fdr(row.get('FDR'))
        
        # Фильтр достоверности (FDR < 0.05)
        if fdr_val is not None and fdr_val < 0.05:
            logfc = row.get('logFC')
            target_id = str(row.get('target_id')).strip() if pd.notna(row.get('target_id')) else ""
            locus = str(row.get('locus')).strip() if pd.notna(row.get('locus')) else ""
            
            if target_id:
                mapping[target_id] = (logfc, fdr_val)
            if locus:
                mapping[locus] = (logfc, fdr_val)
                
    return mapping

def main():
    if len(sys.argv) != 5:
        print("Ошибка: Неверное количество аргументов.")
        print("Использование: python DEGTF.py <таблица TF.xlsx> <таблица DEG1.xlsx> <таблица DEG2.xlsx> <выходная таблица.xlsx>")
        sys.exit(1)

    tf_file = sys.argv[1]
    deg1_file = sys.argv[2]
    deg2_file = sys.argv[3]
    output_file = sys.argv[4]

    for file in [tf_file, deg1_file, deg2_file]:
        if not os.path.exists(file):
            print(f"Ошибка: Файл '{file}' не найден.")
            sys.exit(1)

    print("Загрузка данных...")
    try:
        df_tf = pd.read_excel(tf_file)
        df_deg1 = pd.read_excel(deg1_file)
        df_deg2 = pd.read_excel(deg2_file)
    except Exception as e:
        print(f"Ошибка при чтении файлов Excel: {e}")
        sys.exit(1)

    print("Обработка таблиц DEG и фильтрация по FDR < 0.05...")
    deg1_map = build_deg_mapping(df_deg1)
    deg2_map = build_deg_mapping(df_deg2)

    # Подготовка пустых столбцов
    df_tf['24h'] = None
    df_tf['FDR24'] = None
    df_tf['72h'] = None
    df_tf['FDR72'] = None

    matched_both_count = 0
    total_tf_rows = len(df_tf)

    print("Сопоставление данных (строгое пересечение: DEG1 AND DEG2)...")
    for index, row in df_tf.iterrows():
        locus_tag = str(row.get('Locus Tag')).strip()
        
        if pd.isna(row.get('Locus Tag')) or locus_tag == "nan":
            continue

        # СТРОГОЕ УСЛОВИЕ: Locus Tag должен быть в ОБОИХ словарях
        if locus_tag in deg1_map and locus_tag in deg2_map:
            # Переносим данные из DEG1
            df_tf.at[index, '24h'] = deg1_map[locus_tag][0]
            df_tf.at[index, 'FDR24'] = deg1_map[locus_tag][1]
            # Переносим данные из DEG2
            df_tf.at[index, '72h'] = deg2_map[locus_tag][0]
            df_tf.at[index, 'FDR72'] = deg2_map[locus_tag][1]
            
            matched_both_count += 1

    print("Сохранение результатов...")
    try:
        df_tf.to_excel(output_file, index=False)
        print(f"Файл успешно сохранен как '{output_file}'")
    except Exception as e:
        print(f"Ошибка при сохранении файла: {e}")
        sys.exit(1)

    # Обновленная статистика
    print("\n--- СТАТИСТИКА НАПОЛНЕННОСТИ ---")
    print(f"Всего транскрипционных факторов в таблице TF: {total_tf_rows}")
    print(f"Успешно перенесено (найдено в обеих таблицах с FDR < 0.05): {matched_both_count} шт.")
    print("--------------------------------")

if __name__ == "__main__":
    main()
