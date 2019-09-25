object FormUnsteady: TFormUnsteady
  Left = 0
  Top = 0
  ClientHeight = 392
  ClientWidth = 194
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  OnClose = FormClose
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object RadioGroup1: TRadioGroup
    Left = 0
    Top = 0
    Width = 185
    Height = 49
    Caption = 'Transient Settings'
    Columns = 2
    ItemIndex = 0
    Items.Strings = (
      'steady'
      'unsteady')
    TabOrder = 0
    OnClick = RadioGroup1Click
  end
  object PanelTime: TPanel
    Left = 0
    Top = 111
    Width = 185
    Height = 98
    TabOrder = 1
    Visible = False
    object LabelEndTime: TLabel
      Left = 8
      Top = 20
      Width = 43
      Height = 13
      Caption = 'End Time'
    end
    object LabelTimeUnion: TLabel
      Left = 156
      Top = 22
      Width = 5
      Height = 13
      Caption = 's'
    end
    object EditTime: TEdit
      Left = 57
      Top = 14
      Width = 93
      Height = 21
      TabOrder = 0
      Text = '1'
    end
    object ComboBoxTimeStep: TComboBox
      Left = 5
      Top = 41
      Width = 145
      Height = 21
      ItemIndex = 0
      TabOrder = 1
      Text = 'Linear'
      OnChange = ComboBoxTimeStepChange
      Items.Strings = (
        'Linear'
        'Square Wave'
        'Square Wave APPARAT'
        'hot cold (double linear)')
    end
    object ButtonTimeStepLaw: TButton
      Left = 8
      Top = 68
      Width = 75
      Height = 25
      Caption = 'Edit'
      TabOrder = 2
      Visible = False
      OnClick = ButtonTimeStepLawClick
    end
  end
  object GroupBox1: TGroupBox
    Left = 9
    Top = 215
    Width = 169
    Height = 177
    Caption = 'Additional Setting'
    TabOrder = 2
    object CheckBoxdonttec360: TCheckBox
      Left = 16
      Top = 25
      Width = 137
      Height = 17
      Caption = 'don'#39't start tecplot360'
      TabOrder = 0
    end
    object CheckBoxonlysolidvisible: TCheckBox
      Left = 16
      Top = 48
      Width = 113
      Height = 17
      Caption = 'only solid visible'
      TabOrder = 1
    end
    object CheckBoxreconstruct: TCheckBox
      Left = 16
      Top = 71
      Width = 97
      Height = 17
      Caption = 'reconstruct'
      Checked = True
      State = cbChecked
      TabOrder = 2
    end
    object CheckBoxAnimationFields: TCheckBox
      Left = 16
      Top = 94
      Width = 121
      Height = 17
      Caption = 'save animation fields'
      TabOrder = 3
    end
    object CheckBoxCylinderToPrism: TCheckBox
      Left = 16
      Top = 117
      Width = 121
      Height = 17
      Caption = 'Cylinder_to_Prism'
      TabOrder = 4
    end
    object ButtonCalc: TButton
      Left = 33
      Top = 140
      Width = 120
      Height = 25
      Caption = 'Calculation start'
      TabOrder = 5
      OnClick = ButtonCalcClick
    end
  end
  object GroupBoxNumberIterationsSimpleAlgorithm: TGroupBox
    Left = 1
    Top = 55
    Width = 185
    Height = 50
    Caption = '# iterations SIMPLE algorithm'
    TabOrder = 3
    object ComboBoxNumberIterationsSIMPLEAlgoritm: TComboBox
      Left = 15
      Top = 16
      Width = 106
      Height = 21
      ItemIndex = 0
      TabOrder = 0
      Text = 'default'
      Items.Strings = (
        'default'
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
        '110'
        '111'
        '112'
        '113'
        '114'
        '115'
        '116'
        '117'
        '118'
        '119'
        '120'
        '121'
        '122'
        '123'
        '124'
        '125'
        '126'
        '127'
        '128'
        '129'
        '130'
        '131'
        '132'
        '133'
        '134'
        '135'
        '136'
        '137'
        '138'
        '139'
        '140'
        '141'
        '142'
        '143'
        '144'
        '145'
        '146'
        '147'
        '148'
        '149'
        '150'
        '151'
        '152'
        '153'
        '154'
        '155'
        '156'
        '157'
        '158'
        '159'
        '160'
        '161'
        '162'
        '163'
        '164'
        '165'
        '166'
        '167'
        '168'
        '169'
        '170'
        '171'
        '172'
        '173'
        '174'
        '175'
        '176'
        '177'
        '178'
        '179'
        '180'
        '181'
        '182'
        '183'
        '184'
        '185'
        '186'
        '187'
        '188'
        '189'
        '190'
        '191'
        '192'
        '193'
        '194'
        '195'
        '196'
        '197'
        '198'
        '199'
        '200'
        '201'
        '202'
        '203'
        '204'
        '205'
        '206'
        '207'
        '208'
        '209'
        '210'
        '211'
        '212'
        '213'
        '214'
        '215'
        '216'
        '217'
        '218'
        '219'
        '220'
        '221'
        '222'
        '223'
        '224'
        '225'
        '226'
        '227'
        '228'
        '229'
        '230'
        '231'
        '232'
        '233'
        '234'
        '235'
        '236'
        '237'
        '238'
        '239'
        '240'
        '241'
        '242'
        '243'
        '244'
        '245'
        '246'
        '247'
        '248'
        '249'
        '250'
        '251'
        '252'
        '253'
        '254'
        '255'
        '256'
        '257'
        '258'
        '259'
        '260'
        '261'
        '262'
        '263'
        '264'
        '265'
        '266'
        '267'
        '268'
        '269'
        '270'
        '271'
        '272'
        '273'
        '274'
        '275'
        '276'
        '277'
        '278'
        '279'
        '280'
        '281'
        '282'
        '283'
        '284'
        '285'
        '286'
        '287'
        '288'
        '289'
        '290'
        '291'
        '292'
        '293'
        '294'
        '295'
        '296'
        '297'
        '298'
        '299'
        '300')
    end
  end
end
