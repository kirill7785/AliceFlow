object Form_amg_manager: TForm_amg_manager
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'amg manager (launcher)'
  ClientHeight = 633
  ClientWidth = 494
  Color = clMoneyGreen
  CustomTitleBar.CaptionAlignment = taCenter
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Panel1: TPanel
    Left = 0
    Top = 0
    Width = 494
    Height = 633
    Color = clMoneyGreen
    ParentBackground = False
    TabOrder = 0
    object Label1: TLabel
      Left = 8
      Top = 88
      Width = 132
      Height = 13
      Caption = '2. maximum reduced  levels'
    end
    object Label2: TLabel
      Left = 13
      Top = 223
      Width = 73
      Height = 13
      Caption = '7. interpolation'
    end
    object Label3: TLabel
      Left = 11
      Top = 118
      Width = 54
      Height = 13
      Caption = '3. nFinnest'
    end
    object Label4: TLabel
      Left = 8
      Top = 137
      Width = 121
      Height = 13
      Caption = '4. number pre-smoothing'
    end
    object Label5: TLabel
      Left = 9
      Top = 169
      Width = 120
      Height = 13
      Caption = '5. number post-smothing'
    end
    object Label6: TLabel
      Left = 10
      Top = 196
      Width = 130
      Height = 13
      Caption = '6. memory size n*sizeof(A)'
    end
    object Label7: TLabel
      Left = 12
      Top = 1
      Width = 366
      Height = 13
      Caption = 
        'Settings only for home (original) code RUMBA v0.14 for solve sys' +
        'tem A*x=b'
    end
    object Label8: TLabel
      Left = 10
      Top = 39
      Width = 147
      Height = 13
      Caption = '1. strong connection threshold'
    end
    object Labelmagic: TLabel
      Left = 9
      Top = 319
      Width = 94
      Height = 13
      Caption = '9. F-to-F  threshold'
    end
    object Label9: TLabel
      Left = 8
      Top = 357
      Width = 131
      Height = 13
      Caption = '10. Relaxation (Smoothing)'
    end
    object Label10: TLabel
      Left = 43
      Top = 58
      Width = 60
      Height = 13
      Caption = '(0.23 .. 0.9)'
    end
    object Label11: TLabel
      Left = 92
      Top = 223
      Width = 60
      Height = 13
      Caption = '( 1,2,4 .. 6 )'
    end
    object Label13: TLabel
      Left = 43
      Top = 338
      Width = 60
      Height = 13
      Caption = '(0.35 .. 0.4)'
    end
    object Label14: TLabel
      Left = 164
      Top = 20
      Width = 62
      Height = 13
      Caption = 'Temperature'
    end
    object Label15: TLabel
      Left = 265
      Top = 20
      Width = 30
      Height = 13
      Caption = 'Speed'
    end
    object Label16: TLabel
      Left = 333
      Top = 20
      Width = 42
      Height = 13
      Caption = 'Pressure'
    end
    object Label17: TLabel
      Left = 7
      Top = 408
      Width = 142
      Height = 13
      Caption = '11. amg splitting (coarsening)'
    end
    object Label18: TLabel
      Left = 7
      Top = 435
      Width = 75
      Height = 13
      Caption = '12. stabilization'
    end
    object Label19: TLabel
      Left = 11
      Top = 278
      Width = 138
      Height = 13
      Caption = '8. truncation of interpolation'
    end
    object Label20: TLabel
      Left = 20
      Top = 20
      Width = 46
      Height = 13
      Caption = 'variables '
    end
    object Label21: TLabel
      Left = 8
      Top = 493
      Width = 58
      Height = 13
      Caption = '14. print log'
    end
    object Label26: TLabel
      Left = 9
      Top = 470
      Width = 57
      Height = 13
      Caption = '13. selector'
    end
    object LabelStress: TLabel
      Left = 420
      Top = 20
      Width = 30
      Height = 13
      Caption = 'Stress'
    end
    object Label22: TLabel
      Left = 9
      Top = 580
      Width = 88
      Height = 13
      Caption = '16. Matrix portrait'
    end
    object Buttondefault: TButton
      Left = 68
      Top = 599
      Width = 217
      Height = 25
      Caption = 'return default parameters'
      TabOrder = 0
      OnClick = ButtondefaultClick
    end
    object CheckBoxprintlogTemperature: TCheckBox
      Left = 156
      Top = 492
      Width = 62
      Height = 17
      Caption = 'print'
      Checked = True
      State = cbChecked
      TabOrder = 1
    end
    object CheckBoxprintlogSpeed: TCheckBox
      Left = 246
      Top = 492
      Width = 49
      Height = 17
      Caption = 'print'
      Checked = True
      State = cbChecked
      TabOrder = 2
    end
    object CheckBoxprintlogPressure: TCheckBox
      Left = 326
      Top = 492
      Width = 49
      Height = 17
      Caption = 'print'
      Checked = True
      State = cbChecked
      TabOrder = 3
    end
    object CheckBoxprintlogStress: TCheckBox
      Left = 420
      Top = 492
      Width = 41
      Height = 17
      Caption = 'print'
      Checked = True
      State = cbChecked
      TabOrder = 4
    end
    object GroupBox1: TGroupBox
      Left = 2
      Top = 525
      Width = 492
      Height = 49
      Caption = '15. C-F decomposition Algorithms and Data Structure'
      TabOrder = 5
      object ComboBoxCFalgorithmandDataStruct_Temperature: TComboBox
        Left = 136
        Top = 17
        Width = 91
        Height = 21
        ItemIndex = 2
        TabOrder = 0
        Text = 'Binary Heap [1964]'
        Items.Strings = (
          'AVL Tree [1962]'
          'SPLAY Tree [1985]'
          'Binary Heap [1964]'
          'Treap [1989]'
          'Red-Black Tree [1972]'
          'Fibonacci Heap [1984]'
          'van Emde Boas Tree [1977]')
      end
      object ComboBoxCFalgorithmandDataStruct_Speed: TComboBox
        Left = 233
        Top = 17
        Width = 81
        Height = 21
        ItemIndex = 2
        TabOrder = 1
        Text = 'Binary Heap [1964]'
        Items.Strings = (
          'AVL Tree [1962]'
          'SPLAY Tree [1985]'
          'Binary Heap [1964]'
          'Treap [1989]'
          'Red-Black Tree [1972]'
          'Fibonacci Heap [1984]'
          'van Emde Boas Tree [1977]')
      end
      object ComboBoxCFalgorithmandDataStruct_Pressure: TComboBox
        Left = 320
        Top = 17
        Width = 74
        Height = 21
        ItemIndex = 2
        TabOrder = 2
        Text = 'Binary Heap [1964]'
        Items.Strings = (
          'AVL Tree [1962]'
          'SPLAY Tree [1985]'
          'Binary Heap [1964]'
          'Treap [1989]'
          'Red-Black Tree [1972]'
          'Fibonacci Heap [1984]'
          'van Emde Boas Tree [1977]')
      end
      object ComboBoxCFalgorithmandDataStruct_Stress: TComboBox
        Left = 400
        Top = 17
        Width = 77
        Height = 21
        ItemIndex = 2
        TabOrder = 3
        Text = 'Binary Heap [1964]'
        Items.Strings = (
          'AVL Tree [1962]'
          'SPLAY Tree [1985]'
          'Binary Heap [1964]'
          'Treap [1989]'
          'Red-Black Tree [1972]'
          'Fibonacci Heap [1984]'
          'van Emde Boas Tree [1977]')
      end
    end
    object CheckBoxTemperatureMatrixPortrait: TCheckBox
      Left = 156
      Top = 576
      Width = 50
      Height = 17
      Caption = 'print'
      TabOrder = 6
    end
    object CheckBoxSpeedMatrixPortrait: TCheckBox
      Left = 243
      Top = 576
      Width = 49
      Height = 17
      Caption = 'print'
      TabOrder = 7
    end
    object CheckBoxPressureMatrixPortrait: TCheckBox
      Left = 331
      Top = 576
      Width = 50
      Height = 17
      Caption = 'print'
      TabOrder = 8
    end
    object CheckBoxStressMatrixPortrait: TCheckBox
      Left = 420
      Top = 576
      Width = 49
      Height = 17
      Caption = 'print'
      TabOrder = 9
    end
    object ComboBoxSort: TComboBox
      Left = 332
      Top = 250
      Width = 145
      Height = 21
      ItemIndex = 0
      TabOrder = 10
      Text = 'COUNTING SORT'
      Items.Strings = (
        'COUNTING SORT'
        'QUICK SORT'
        'HEAP SORT'
        'TIM SORT'
        'in development')
    end
    object CheckBoxDiagonalDominant: TCheckBox
      Left = 162
      Top = 254
      Width = 111
      Height = 17
      Caption = 'diagonal dominance'
      Checked = True
      State = cbChecked
      TabOrder = 11
    end
    object CheckBoxStrongTranspose: TCheckBox
      Left = 156
      Top = 469
      Width = 117
      Height = 17
      Caption = 'Strong Transpose'
      Checked = True
      State = cbChecked
      TabOrder = 12
    end
    object PanelTemperature1: TPanel
      Left = 158
      Top = 39
      Width = 81
      Height = 209
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 13
      object Editthreshold: TEdit
        Left = 8
        Top = 20
        Width = 66
        Height = 21
        TabOrder = 0
        Text = '0.5'
      end
      object ComboBoxmaximumreducedlevels: TComboBox
        Left = 8
        Top = 47
        Width = 66
        Height = 21
        ItemIndex = 0
        TabOrder = 1
        Text = '0'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6'
          '7'
          '8'
          '9'
          '10'
          '11'
          '12'
          '13'
          '14'
          '15'
          '16'
          '17'
          '18'
          '19'
          '20')
      end
      object ComboBoxnFinnest: TComboBox
        Left = 8
        Top = 74
        Width = 66
        Height = 21
        ItemIndex = 1
        TabOrder = 2
        Text = '2'
        Items.Strings = (
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object ComboBoxnumberpresmothers: TComboBox
        Left = 8
        Top = 101
        Width = 66
        Height = 21
        ItemIndex = 1
        TabOrder = 3
        Text = '1'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object ComboBoxnumberpostsweeps: TComboBox
        Left = 8
        Top = 128
        Width = 66
        Height = 21
        ItemIndex = 2
        TabOrder = 4
        Text = '2'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object ComboBoxmemorysize: TComboBox
        Left = 8
        Top = 155
        Width = 66
        Height = 21
        ItemIndex = 2
        TabOrder = 5
        Text = '4'
        Items.Strings = (
          '2'
          '3'
          '4'
          '5'
          '6'
          '7'
          '8'
          '9'
          '10'
          '11'
          '12'
          '13'
          '14'
          '15'
          '16'
          '17'
          '18'
          '19'
          '20'
          '21'
          '22'
          '23'
          '24'
          '25'
          '26'
          '27'
          '28'
          '29'
          '30'
          '31'
          '32'
          '33'
          '34'
          '35'
          '36'
          '37'
          '38'
          '39'
          '40'
          '41'
          '42'
          '43'
          '44'
          '45'
          '46'
          '47'
          '48'
          '49'
          '50'
          '51'
          '52'
          '53'
          '54'
          '55'
          '56'
          '57'
          '58'
          '59'
          '60'
          '61'
          '62'
          '63'
          '64'
          '65'
          '66'
          '67'
          '68'
          '69'
          '70'
          '71'
          '72'
          '73'
          '74'
          '75'
          '76'
          '77'
          '78'
          '79'
          '80'
          '81'
          '82'
          '83'
          '84'
          '85'
          '86'
          '87'
          '88'
          '89'
          '90'
          '91'
          '92'
          '93'
          '94'
          '95'
          '96'
          '97'
          '98'
          '99'
          '100')
      end
      object ComboBoxinterpolation: TComboBox
        Left = 8
        Top = 182
        Width = 66
        Height = 21
        ItemIndex = 1
        TabOrder = 6
        Text = '2. J.W.RUGE and K.STUBEN [1987] p.30 (102) AMG1R5'
        Items.Strings = (
          '1. JACOBI distance 2 interpolation'
          '2. J.W.RUGE and K.STUBEN [1987] p.30 (102) AMG1R5'
          '3  Longe Range interpolation'
          '4.  for a long time usage (lite 1)'
          '5.  lite version 4'
          '6.  square lite version 4')
      end
      object CheckBoxbautoThresholdTemperature: TCheckBox
        Left = 5
        Top = 0
        Width = 97
        Height = 17
        Caption = 'auto'
        TabOrder = 7
        OnClick = CheckBoxbautoThresholdTemperatureClick
      end
    end
    object PanelSpeed1: TPanel
      Left = 238
      Top = 39
      Width = 89
      Height = 209
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 14
      object EditthresholdSpeed: TEdit
        Left = 9
        Top = 23
        Width = 72
        Height = 21
        TabOrder = 0
        Text = '0.5'
      end
      object ComboBoxmaximumreducedlevelsSpeed: TComboBox
        Left = 10
        Top = 46
        Width = 71
        Height = 21
        ItemIndex = 0
        TabOrder = 1
        Text = '0'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6'
          '7'
          '8'
          '9'
          '10'
          '11'
          '12'
          '13'
          '14'
          '15'
          '16'
          '17'
          '18'
          '19'
          '20')
      end
      object ComboBoxnFinnestSpeed: TComboBox
        Left = 10
        Top = 73
        Width = 71
        Height = 21
        ItemIndex = 1
        TabOrder = 2
        Text = '2'
        Items.Strings = (
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object ComboBoxnumberpresmothersSpeed: TComboBox
        Left = 10
        Top = 100
        Width = 71
        Height = 21
        ItemIndex = 1
        TabOrder = 3
        Text = '1'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object ComboBoxnumberpostsweepsSpeed: TComboBox
        Left = 10
        Top = 127
        Width = 71
        Height = 21
        ItemIndex = 2
        TabOrder = 4
        Text = '2'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object ComboBoxmemorysizeSpeed: TComboBox
        Left = 10
        Top = 154
        Width = 71
        Height = 21
        ItemIndex = 2
        TabOrder = 5
        Text = '4'
        Items.Strings = (
          '2'
          '3'
          '4'
          '5'
          '6'
          '7'
          '8'
          '9'
          '10'
          '11'
          '12'
          '13'
          '14'
          '15'
          '16'
          '17'
          '18'
          '19'
          '20'
          '21'
          '22'
          '23'
          '24'
          '25'
          '26'
          '27'
          '28'
          '29'
          '30'
          '31'
          '32'
          '33'
          '34'
          '35'
          '36'
          '37'
          '38'
          '39'
          '40'
          '41'
          '42'
          '43'
          '44'
          '45'
          '46'
          '47'
          '48'
          '49'
          '50'
          '51'
          '52'
          '53'
          '54'
          '55'
          '56'
          '57'
          '58'
          '59'
          '60'
          '61'
          '62'
          '63'
          '64'
          '65'
          '66'
          '67'
          '68'
          '69'
          '70'
          '71'
          '72'
          '73'
          '74'
          '75'
          '76'
          '77'
          '78'
          '79'
          '80'
          '81'
          '82'
          '83'
          '84'
          '85'
          '86'
          '87'
          '88'
          '89'
          '90'
          '91'
          '92'
          '93'
          '94'
          '95'
          '96'
          '97'
          '98'
          '99'
          '100')
      end
      object ComboBoxInterpolationSpeed: TComboBox
        Left = 10
        Top = 181
        Width = 71
        Height = 21
        ItemIndex = 1
        TabOrder = 6
        Text = '2. J.W.RUGE and K.STUBEN [1987] p.30 (102) AMG1R5'
        Items.Strings = (
          '1. JACOBI distance 2 interpolation'
          '2. J.W.RUGE and K.STUBEN [1987] p.30 (102) AMG1R5'
          '3.  Longe Range interpolation'
          '4.  for a long time usage (lite 1)'
          '5.  lite version 4'
          '6.  square lite version 4')
      end
      object CheckBoxbautoThresholdSpeed: TCheckBox
        Left = 16
        Top = 0
        Width = 97
        Height = 17
        Caption = 'auto'
        TabOrder = 7
        OnClick = CheckBoxbautoThresholdSpeedClick
      end
    end
    object PanelPressure1: TPanel
      Left = 325
      Top = 39
      Width = 81
      Height = 209
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 15
      object EditthresholdPressure: TEdit
        Left = 8
        Top = 22
        Width = 65
        Height = 21
        TabOrder = 0
        Text = '0.5'
      end
      object ComboBoxmaximumreducedlevelsPressure: TComboBox
        Left = 8
        Top = 49
        Width = 65
        Height = 21
        ItemIndex = 0
        TabOrder = 1
        Text = '0'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6'
          '7'
          '8'
          '9'
          '10'
          '11'
          '12'
          '13'
          '14'
          '15'
          '16'
          '17'
          '18'
          '19'
          '20')
      end
      object ComboBoxnFinnestPressure: TComboBox
        Left = 8
        Top = 76
        Width = 65
        Height = 21
        ItemIndex = 1
        TabOrder = 2
        Text = '2'
        Items.Strings = (
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object ComboBoxnumberpresmothersPressure: TComboBox
        Left = 8
        Top = 103
        Width = 65
        Height = 21
        ItemIndex = 1
        TabOrder = 3
        Text = '1'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object ComboBoxnumberpostsweepsPressure: TComboBox
        Left = 7
        Top = 130
        Width = 66
        Height = 21
        ItemIndex = 2
        TabOrder = 4
        Text = '2'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object ComboBoxmemorysizePressure: TComboBox
        Left = 7
        Top = 157
        Width = 66
        Height = 21
        ItemIndex = 2
        TabOrder = 5
        Text = '4'
        Items.Strings = (
          '2'
          '3'
          '4'
          '5'
          '6'
          '7'
          '8'
          '9'
          '10'
          '11'
          '12'
          '13'
          '14'
          '15'
          '16'
          '17'
          '18'
          '19'
          '20'
          '21'
          '22'
          '23'
          '24'
          '25'
          '26'
          '27'
          '28'
          '29'
          '30'
          '31'
          '32'
          '33'
          '34'
          '35'
          '36'
          '37'
          '38'
          '39'
          '40'
          '41'
          '42'
          '43'
          '44'
          '45'
          '46'
          '47'
          '48'
          '49'
          '50'
          '51'
          '52'
          '53'
          '54'
          '55'
          '56'
          '57'
          '58'
          '59'
          '60'
          '61'
          '62'
          '63'
          '64'
          '65'
          '66'
          '67'
          '68'
          '69'
          '70'
          '71'
          '72'
          '73'
          '74'
          '75'
          '76'
          '77'
          '78'
          '79'
          '80'
          '81'
          '82'
          '83'
          '84'
          '85'
          '86'
          '87'
          '88'
          '89'
          '90'
          '91'
          '92'
          '93'
          '94'
          '95'
          '96'
          '97'
          '98'
          '99'
          '100')
      end
      object ComboBoxinterpolationPressure: TComboBox
        Left = 7
        Top = 184
        Width = 66
        Height = 21
        ItemIndex = 1
        TabOrder = 6
        Text = '2. J.W.RUGE and K.STUBEN [1987] p.30 (102) AMG1R5'
        Items.Strings = (
          '1. JACOBI distance 2 interpolation'
          '2. J.W.RUGE and K.STUBEN [1987] p.30 (102) AMG1R5'
          '3.  Longe Range interpolation'
          '4.  for a long time usage (lite 1)'
          '5.  lite version 4'
          '6.  square lite version 4')
      end
      object CheckBoxbthresholdPressureAuto: TCheckBox
        Left = 16
        Top = 0
        Width = 97
        Height = 17
        Caption = 'auto'
        TabOrder = 7
        OnClick = CheckBoxbthresholdPressureAutoClick
      end
    end
    object PanelStress1: TPanel
      Left = 404
      Top = 40
      Width = 81
      Height = 208
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 16
      object EditthresholdStress: TEdit
        Left = 8
        Top = 21
        Width = 65
        Height = 21
        TabOrder = 0
        Text = '0.5'
      end
      object ComboBoxmaximumreducedlevelsStress: TComboBox
        Left = 8
        Top = 48
        Width = 65
        Height = 21
        ItemIndex = 0
        TabOrder = 1
        Text = '0'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6'
          '7'
          '8'
          '9'
          '10'
          '11'
          '12'
          '13'
          '14'
          '15'
          '16'
          '17'
          '18'
          '19'
          '20')
      end
      object ComboBoxnFinnestStress: TComboBox
        Left = 8
        Top = 75
        Width = 65
        Height = 21
        ItemIndex = 1
        TabOrder = 2
        Text = '2'
        Items.Strings = (
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object ComboBoxnumberpresmoothersStress: TComboBox
        Left = 8
        Top = 102
        Width = 65
        Height = 21
        ItemIndex = 1
        TabOrder = 3
        Text = '1'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object ComboBoxnumberpostsweepsStress: TComboBox
        Left = 8
        Top = 129
        Width = 65
        Height = 21
        ItemIndex = 2
        TabOrder = 4
        Text = '2'
        Items.Strings = (
          '0'
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object ComboBoxmemorysizeStress: TComboBox
        Left = 8
        Top = 156
        Width = 65
        Height = 21
        ItemIndex = 5
        TabOrder = 5
        Text = '7'
        Items.Strings = (
          '2'
          '3'
          '4'
          '5'
          '6'
          '7'
          '8'
          '9'
          '10'
          '11'
          '12'
          '13'
          '14'
          '15'
          '16'
          '17'
          '18'
          '19'
          '20'
          '21'
          '22'
          '23'
          '24'
          '25'
          '26'
          '27'
          '28'
          '29'
          '30'
          '31'
          '32'
          '33'
          '34'
          '35'
          '36'
          '37'
          '38'
          '39'
          '40'
          '41'
          '42'
          '43'
          '44'
          '45'
          '46'
          '47'
          '48'
          '49'
          '50'
          '51'
          '52'
          '53'
          '54'
          '55'
          '56'
          '57'
          '58'
          '59'
          '60'
          '61'
          '62'
          '63'
          '64'
          '65'
          '66'
          '67'
          '68'
          '69'
          '70'
          '71'
          '72'
          '73'
          '74'
          '75'
          '76'
          '77'
          '78'
          '79'
          '80'
          '81'
          '82'
          '83'
          '84'
          '85'
          '86'
          '87'
          '88'
          '89'
          '90'
          '91'
          '92'
          '93'
          '94'
          '95'
          '96'
          '97'
          '98'
          '99'
          '100'
          '101'
          '102'
          '103'
          '104'
          '105'
          '106'
          '107'
          '108'
          '109'
          '110')
      end
      object ComboBoxinterpollationStress: TComboBox
        Left = 8
        Top = 183
        Width = 65
        Height = 21
        ItemIndex = 1
        TabOrder = 6
        Text = '2. J.W.RUGE and K.STUBEN [1987] p.30 (102) AMG1R5'
        Items.Strings = (
          '1. JACOBI distance 2 interpolation'
          '2. J.W.RUGE and K.STUBEN [1987] p.30 (102) AMG1R5'
          '3.  Longe Range interpolation'
          '4.  for a long time usage (lite 1)'
          '5.  lite version 4'
          '6.  square lite version 4')
      end
      object CheckBoxbthresholdautoStress: TCheckBox
        Left = 8
        Top = 0
        Width = 97
        Height = 17
        Caption = 'auto'
        TabOrder = 7
        OnClick = CheckBoxbthresholdautoStressClick
      end
    end
    object PanelTemperature2: TPanel
      Left = 155
      Top = 277
      Width = 81
      Height = 186
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 17
      object CheckBoxtruncationT: TCheckBox
        Left = 8
        Top = 0
        Width = 41
        Height = 17
        Caption = 'off'
        TabOrder = 0
        OnClick = CheckBoxtruncationTClick
      end
      object Edit_truncation_T: TEdit
        Left = 7
        Top = 23
        Width = 66
        Height = 21
        TabOrder = 1
        Text = '0.3'
      end
      object EditmagicT: TEdit
        Left = 8
        Top = 50
        Width = 65
        Height = 21
        TabOrder = 2
        Text = '0.4'
      end
      object ComboBoxsmoothertypeTemperature: TComboBox
        Left = 8
        Top = 77
        Width = 65
        Height = 21
        ItemIndex = 0
        TabOrder = 3
        Text = 'Gauss-Seidel'
        Items.Strings = (
          'Gauss-Seidel'
          'iluk(k== lfil)'
          'Runge-Kutta order = 3'
          'Runge-Kutta order = 5'
          'damped Jacoby'
          'Rouch sor'
          'GMRES'
          'spai-0'
          'Chebyshev P.L.')
      end
      object ComboBoxcoarseningTemp: TComboBox
        Left = 11
        Top = 128
        Width = 64
        Height = 21
        ItemIndex = 10
        TabOrder = 4
        Text = 'PMIS2'
        OnChange = ComboBoxcoarseningTempChange
        Items.Strings = (
          'classical all con'
          'RS 2 all con'
          'classical ST all con'
          'RS 2 ST all con'
          'classical neg con'
          'RS 2 neg con'
          'classical ST neg con'
          'RS 2 ST neg con'
          'PMIS'
          'HMIS'
          'PMIS2'
          'PMIS2_full')
      end
      object ComboBoxStabilizationTemp: TComboBox
        Left = 9
        Top = 155
        Width = 64
        Height = 21
        ItemIndex = 2
        TabOrder = 5
        Text = 'FGMRes [1986]'
        Items.Strings = (
          'none'
          'BiCGStab [1992]'
          'FGMRes [1986]'
          'for_NonLinear_problem')
      end
      object CheckBoxCFReorderTemperature: TCheckBox
        Left = 5
        Top = 105
        Width = 97
        Height = 17
        Caption = 'CF reorder'
        TabOrder = 6
      end
    end
    object PanelSpeed2: TPanel
      Left = 234
      Top = 277
      Width = 81
      Height = 186
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 18
      object CheckBoxtruncationSpeed: TCheckBox
        Left = 6
        Top = 0
        Width = 43
        Height = 17
        Caption = 'off'
        Checked = True
        State = cbChecked
        TabOrder = 0
        OnClick = CheckBoxtruncationSpeedClick
      end
      object Edit_truncation_Speed: TEdit
        Left = 7
        Top = 23
        Width = 66
        Height = 21
        TabOrder = 1
        Text = '0.3'
        Visible = False
      end
      object EditmagicSpeed: TEdit
        Left = 7
        Top = 50
        Width = 66
        Height = 21
        TabOrder = 2
        Text = '0.4'
      end
      object ComboBoxsmoothertypeSpeed: TComboBox
        Left = 7
        Top = 77
        Width = 66
        Height = 21
        ItemIndex = 0
        TabOrder = 3
        Text = 'Gauss-Seidel'
        Items.Strings = (
          'Gauss-Seidel'
          'iluk(k== lfil)'
          'Runge-Kutta order = 3'
          'Runge-Kutta order = 5'
          'damped Jacoby'
          'Rouch sor'
          'GMRES'
          'spai-0')
      end
      object ComboBoxcoarseningSpeed: TComboBox
        Left = 8
        Top = 128
        Width = 66
        Height = 21
        ItemIndex = 10
        TabOrder = 4
        Text = 'PMIS2'
        OnChange = ComboBoxcoarseningSpeedChange
        Items.Strings = (
          'classical all con'
          'RS 2 all con'
          'classical ST all con'
          'RS 2 ST all con'
          'classical neg con'
          'RS 2 neg con'
          'classical ST neg con'
          'RS 2 ST neg con'
          'PMIS'
          'HMIS'
          'PMIS2'
          'PMIS2_full')
      end
      object ComboBoxStabilizationSpeed: TComboBox
        Left = 8
        Top = 155
        Width = 66
        Height = 21
        ItemIndex = 2
        TabOrder = 5
        Text = 'FGMRes [1986]'
        Items.Strings = (
          'none'
          'BiCGStab [1992]'
          'FGMRes [1986]')
      end
      object CheckBoxCFReorderSpeed: TCheckBox
        Left = 8
        Top = 105
        Width = 97
        Height = 17
        Caption = 'CF reorder'
        TabOrder = 6
      end
    end
    object PanelPressure2: TPanel
      Left = 315
      Top = 277
      Width = 81
      Height = 186
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 19
      object CheckBoxtruncationPressure: TCheckBox
        Left = 8
        Top = 0
        Width = 49
        Height = 17
        Caption = 'off'
        Checked = True
        State = cbChecked
        TabOrder = 0
        OnClick = CheckBoxtruncationPressureClick
      end
      object Edit_truncation_Pressure: TEdit
        Left = 7
        Top = 23
        Width = 66
        Height = 21
        TabOrder = 1
        Text = '0.3'
        Visible = False
      end
      object EditmagicPressure: TEdit
        Left = 7
        Top = 50
        Width = 66
        Height = 21
        TabOrder = 2
        Text = '0.4'
      end
      object ComboBoxsmoothertypePressure: TComboBox
        Left = 9
        Top = 77
        Width = 64
        Height = 21
        ItemIndex = 5
        TabOrder = 3
        Text = 'Rouch sor'
        Items.Strings = (
          'Gauss-Seidel'
          'iluk(k== lfil)'
          'Runge-Kutta order=3'
          'Runge-Kutta order=5'
          'damped Jacoby'
          'Rouch sor'
          'GMRES'
          'spai0'
          'Chebyshev P.L.')
      end
      object ComboBoxcoarseningPressure: TComboBox
        Left = 6
        Top = 129
        Width = 65
        Height = 21
        ItemIndex = 10
        TabOrder = 4
        Text = 'PMIS2'
        OnChange = ComboBoxcoarseningPressureChange
        Items.Strings = (
          'classical all con'
          'RS 2 all con'
          'classical ST all con'
          'RS 2 ST all con'
          'classical neg con'
          'RS 2 neg con'
          'classical ST neg con'
          'RS 2 ST neg con'
          'PMIS'
          'HMIS'
          'PMIS2'
          'PMIS2_full')
      end
      object ComboBoxStabilizationPressure: TComboBox
        Left = 7
        Top = 155
        Width = 66
        Height = 21
        ItemIndex = 2
        TabOrder = 5
        Text = 'FGMRes [1986]'
        Items.Strings = (
          'none'
          'BiCGStab [1992]'
          'FGMRes [1986]'
          'CG [1952]')
      end
      object CheckBoxCFReorderPressure: TCheckBox
        Left = 8
        Top = 105
        Width = 97
        Height = 17
        Caption = 'CF reorder'
        TabOrder = 6
      end
    end
    object PanelStress2: TPanel
      Left = 394
      Top = 277
      Width = 81
      Height = 186
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 20
      object CheckBoxtruncationStress: TCheckBox
        Left = 3
        Top = 0
        Width = 30
        Height = 17
        Caption = 'off'
        Checked = True
        State = cbChecked
        TabOrder = 0
        OnClick = CheckBoxtruncationStressClick
      end
      object Edittruncation_Stress: TEdit
        Left = 7
        Top = 213
        Width = 66
        Height = 21
        TabOrder = 1
        Text = '0.2'
        Visible = False
        OnClick = Edittruncation_StressClick
      end
      object EditmagicStress: TEdit
        Left = 8
        Top = 50
        Width = 65
        Height = 21
        TabOrder = 2
        Text = '0.4'
      end
      object ComboBoxsmoothertypeStress: TComboBox
        Left = 4
        Top = 77
        Width = 69
        Height = 21
        ItemIndex = 8
        TabOrder = 3
        Text = 'Chebyshev P.L.'
        Items.Strings = (
          'Gauss-Seidel'
          'iluk(k== lfil)'
          'Runge-Kutta order=3'
          'Runge-Kutta order=5'
          'damped Jacoby'
          'Rouch sor'
          'GMRES'
          'spai-0'
          'Chebyshev P.L.')
      end
      object ComboBoxcoarseningStress: TComboBox
        Left = 8
        Top = 128
        Width = 69
        Height = 21
        ItemIndex = 10
        TabOrder = 4
        Text = 'PMIS2'
        OnChange = ComboBoxcoarseningStressChange
        Items.Strings = (
          'classical all con'
          'RS 2 all con'
          'classical ST all con'
          'RS 2 ST all con'
          'classical neg con'
          'RS 2 neg con'
          'classical ST neg con'
          'RS 2 ST neg con'
          'PMIS'
          'HMIS'
          'PMIS2'
          'PMIS2_full')
      end
      object ComboBoxstabilizationStress: TComboBox
        Left = 8
        Top = 155
        Width = 69
        Height = 21
        ItemIndex = 2
        TabOrder = 5
        Text = 'FGMRes [1986]'
        Items.Strings = (
          'none'
          'BiCGStab [1992]'
          'FGMRes [1986]'
          'CG [1952]')
      end
      object CheckBoxCFReorderStress: TCheckBox
        Left = 8
        Top = 105
        Width = 97
        Height = 17
        Caption = 'CF reorder'
        TabOrder = 6
      end
      object Edit_truncation_Stress: TEdit
        Left = 8
        Top = 23
        Width = 65
        Height = 21
        TabOrder = 7
        Text = '0.3'
      end
    end
    object ButtonRumbaApply: TButton
      Left = 360
      Top = 599
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 21
      OnClick = ButtonRumbaApplyClick
    end
  end
  object ApplicationEvents1: TApplicationEvents
    OnMessage = ApplicationEvents1Message
    Left = 97
    Top = 480
  end
end
