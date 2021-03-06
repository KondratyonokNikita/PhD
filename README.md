<div align="center">
  <H1>
    Свойства теоретико-числовых и криптографических алгоритмов в дедекиндовых кольцах
  </H1>
  Диссертация на соискание ученой степени кандидата физико-математических наук<br><br>
  по специальности <a href="https://vak.gov.by/node/1274">01.01.06</a> – ''Математическая логика, алгебра и теория чисел''<br><br>
  Кондратёнок Никита Васильевич
</div><br>
<div align="center">
  Научный руководитель: д.ф.-м.н. Васьковский Максим Михайлович
</div>

## Планируемое содержание

- Перечень сокращений и (или) условных обозначений
- Введение
- Общая характеристика работы
- <details><summary>Глава 1. <b>Предварительные сведения</b></summary>

  - Определения идеала, простого идеала, максимального идеала, дедекиндова кольца
  - Функция Эйлера в дедекиндовом кольце и ее свойства
  - Теорема Копперсмита
  - Определения нормы, дробной и целой частей, цепочки делений
  - Примеры нормы, пример кольца, где нет цепочки делений с выбором минимального по норме остатка
  - Определение регулярной тройки и формулировка теоремы Кронекера-Валена

</details>

- <details><summary>Глава 2. <b>Алгебраические операции в дедекиндовых кольцах</b> (Тестирование идеалов на простоту в дедекиндовых кольцах)</summary>

  - Факторизация идеалов
    - Использование теоремы Дедекинда для сведения задачи факторизации к целым числам
    - Привести результаты Kofi_Intrinsic factorization of ideals in dedekind domains, где используется вычисление радикала
  - Тестирование идеалов на простоту
    - Аналог критерия Миллера и оценки вероятности успеха
    - Аналог критерия Эйлера и оценки вероятности успеха
    - Детерминированное тестирование на простоту
  - Вычислительная сложность алгебраических операций
    - Вычислительная сложность элементарных операций
    - Сложность вероятностного тестирования на простоту
    - Сложность алгоритма факторизации

</details>

- <details><summary>Глава 3. <b>Теорема Кронекера-Валена в факториальных кольцах</b></summary>

  - Предварительные сведения
  - Теорема Кронекера-Валена в специальном классе факториальных колец
  - Метод проверки принадлежности кольца классу T
    - Определение класса S
    - Доказательство, что S подмножество T
    - Метод проверки принадлежности классу S
    - Примеры из класса S, из T и не из S, не из T
    - Метод проверки принадлежности классу T
  - Теорема Ламе в факториальных кольцах
  - Теорема Кронекера-Валена в кольцах целых алгебраических чисел
    - Определения
    - Алгоритм вычисления наименьшего по норме остатка
    - Вычислительная сложность алгоритма
    - Метод доказательства невыполнимости теоремы Кронекера-Валена
    - Теорема, что для действительных квадратичных норменно-евклидовых колец теорема Кронекера-Валена не выполнена
    - Теорема для всех квадратичных норменно-евклидовых колец.

</details>

- <details><summary>Глава 4. <b>Аналог RSA-криптосистемы в дедекиндовых кольцах</b></summary>

  - Формулировка аналога RSA-криптосистемы
    - Доказательство работоспособности.
    - Ограничения для вычислимости алгоритма (из статьи Petukhova, Tronin_RSA cryptosystem for Dedekind rings)
  - Анализ аналога RSA-криптосистемы
    - Теорема, что если d известно, то N можно разложить с вероятностью не менее 1/2 за лог время. (кажется только для факториальных, так как надо искать НОД(b-1, N))
    - Теорема Винера, что если d маленькое, то его можно вычислить. (для дедекиндовых колец)
    - Метод повторного шифрования. (для дедекиндовых колец)
    - Теорема, что если у нормы p и q одинаковая битовая длина, то их эти нормы можно вычислить. (для дедекиндовых)
    - Теорема, что нельзя иметь одинаковые RSA-модули. (для евклидовых колец)
  - Пример работы криптосистемы в координатных кольцах

</details>

- Заключение
- Библиографический список

## Другая информация

Должно быть около 80 страниц в диссертации. Около 80-100 источников в списке литературы. Весь текст лучше написать заново, а не брать части из диплома или магистерской диссертации.

Приоритет на аккуратность доказательств и формулировок. Важно чтобы не было многозначности.

Прошлая специальность: [05.13.19](https://vak.gov.by/node/1467) – ''Методы и системы защиты информации, информационная безопасность''.

### Структура директорий

- `additional` - дополнительные файлы
- `sources` - исходный код диссертации

### Полезные ссылки

- https://github.com/belgraviton/thesisby
- https://github.com/andriygav/PhDThesis
