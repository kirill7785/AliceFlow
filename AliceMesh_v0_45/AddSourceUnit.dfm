object AddSourceForm: TAddSourceForm
  Left = 157
  Top = 154
  AutoSize = True
  Caption = 'edit properties of the source'
  ClientHeight = 369
  ClientWidth = 333
  Color = clMoneyGreen
  CustomTitleBar.CaptionAlignment = taCenter
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Panelglobalcontainer: TPanel
    Left = 0
    Top = 0
    Width = 333
    Height = 369
    Color = clMoneyGreen
    ParentBackground = False
    TabOrder = 0
    object RadioGroup1: TRadioGroup
      Left = 8
      Top = 8
      Width = 305
      Height = 49
      Caption = 'size and properties source'
      Columns = 4
      ItemIndex = 0
      Items.Strings = (
        'Info'
        'Geometry'
        'Properties'
        'Network')
      TabOrder = 0
      OnClick = RadioGroup1Click
    end
    object Bapply: TButton
      Left = 250
      Top = 336
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 1
      OnClick = BapplyClick
    end
    object Panelinfo: TPanel
      Left = 8
      Top = 83
      Width = 305
      Height = 247
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 2
      object Lname: TLabel
        Left = 16
        Top = 16
        Width = 26
        Height = 13
        Caption = 'name'
      end
      object Ename: TEdit
        Left = 56
        Top = 8
        Width = 105
        Height = 21
        TabOrder = 0
      end
    end
    object PanelGeometry: TPanel
      Left = 8
      Top = 83
      Width = 305
      Height = 247
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 3
      Visible = False
      object GroupBox1: TGroupBox
        Left = 8
        Top = 69
        Width = 257
        Height = 129
        Caption = 'size object'
        Color = clMoneyGreen
        ParentColor = False
        TabOrder = 0
        object LxS: TLabel
          Left = 16
          Top = 32
          Width = 12
          Height = 13
          Caption = 'xS'
        end
        object LyS: TLabel
          Left = 16
          Top = 64
          Width = 12
          Height = 13
          Caption = 'yS'
        end
        object LzS: TLabel
          Left = 16
          Top = 96
          Width = 12
          Height = 13
          Caption = 'zS'
        end
        object LxE: TLabel
          Left = 128
          Top = 32
          Width = 12
          Height = 13
          Caption = 'xE'
        end
        object LyE: TLabel
          Left = 128
          Top = 64
          Width = 12
          Height = 13
          Caption = 'yE'
        end
        object LzE: TLabel
          Left = 128
          Top = 96
          Width = 12
          Height = 13
          Caption = 'zE'
          Visible = False
        end
        object ExS: TEdit
          Left = 40
          Top = 24
          Width = 73
          Height = 21
          TabOrder = 0
        end
        object EyS: TEdit
          Left = 40
          Top = 56
          Width = 73
          Height = 21
          TabOrder = 1
        end
        object EzS: TEdit
          Left = 40
          Top = 88
          Width = 73
          Height = 21
          TabOrder = 2
        end
        object ExE: TEdit
          Left = 152
          Top = 24
          Width = 81
          Height = 21
          TabOrder = 3
        end
        object EyE: TEdit
          Left = 152
          Top = 56
          Width = 81
          Height = 21
          TabOrder = 4
        end
        object EzE: TEdit
          Left = 152
          Top = 88
          Width = 81
          Height = 21
          TabOrder = 5
          Visible = False
        end
      end
      object RadioGroupPlane: TRadioGroup
        Left = 8
        Top = 8
        Width = 225
        Height = 49
        Caption = 'plane'
        Color = clMoneyGreen
        Columns = 3
        ItemIndex = 0
        Items.Strings = (
          'X-Y'
          'X-Z'
          'Y-Z')
        ParentColor = False
        TabOrder = 1
        OnClick = RadioGroupPlaneClick
      end
    end
    object PanelProperties: TPanel
      Left = 0
      Top = 83
      Width = 305
      Height = 247
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 4
      Visible = False
      object GBpowerdef: TGroupBox
        Left = 0
        Top = 8
        Width = 305
        Height = 157
        Caption = 'power define'
        Color = clMoneyGreen
        ParentColor = False
        TabOrder = 0
        object Label1: TLabel
          Left = 16
          Top = 128
          Width = 29
          Height = 13
          Caption = 'power'
        end
        object LW: TLabel
          Left = 126
          Top = 122
          Width = 11
          Height = 13
          Caption = 'W'
        end
        object RGpowertype: TRadioGroup
          Left = 8
          Top = 24
          Width = 145
          Height = 89
          Caption = 'power type'
          ItemIndex = 0
          Items.Strings = (
            'const'
            'temperature depend')
          TabOrder = 0
          OnClick = RGpowertypeClick
        end
        object Ptempdefloc: TPanel
          Left = 159
          Top = 28
          Width = 137
          Height = 97
          Color = clMoneyGreen
          TabOrder = 1
          object Label2: TLabel
            Left = 8
            Top = 16
            Width = 34
            Height = 13
            Caption = 'id table'
          end
          object Label3: TLabel
            Left = 8
            Top = 48
            Width = 99
            Height = 13
            Caption = 'operating offset drain'
          end
          object CBtableid: TComboBox
            Left = 48
            Top = 8
            Width = 73
            Height = 21
            TabOrder = 0
          end
          object EOperoffsetdrain: TEdit
            Left = 8
            Top = 64
            Width = 105
            Height = 21
            TabOrder = 1
          end
        end
        object Epower: TEdit
          Left = 51
          Top = 119
          Width = 65
          Height = 21
          TabOrder = 2
        end
      end
    end
    object PanelNetwork: TPanel
      Left = 0
      Top = 80
      Width = 313
      Height = 250
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 5
      Visible = False
      object GroupBoxinx: TGroupBox
        Left = 16
        Top = 17
        Width = 185
        Height = 49
        Caption = 'inx'
        TabOrder = 0
        object ComboBoxinx: TComboBox
          Left = 24
          Top = 11
          Width = 145
          Height = 21
          ItemIndex = 0
          TabOrder = 0
          Text = '1'
          Items.Strings = (
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
            '100')
        end
      end
      object GroupBoxiny: TGroupBox
        Left = 16
        Top = 74
        Width = 185
        Height = 50
        Caption = 'iny'
        TabOrder = 1
        object ComboBoxiny: TComboBox
          Left = 24
          Top = 16
          Width = 145
          Height = 21
          ItemIndex = 0
          TabOrder = 0
          Text = '1'
          Items.Strings = (
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
            '100')
        end
      end
      object GroupBoxinz: TGroupBox
        Left = 16
        Top = 130
        Width = 185
        Height = 55
        Caption = 'inz'
        TabOrder = 2
        object ComboBoxinz: TComboBox
          Left = 24
          Top = 22
          Width = 145
          Height = 21
          ItemIndex = 0
          TabOrder = 0
          Text = '1'
          Items.Strings = (
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
            '100')
        end
      end
    end
  end
end
